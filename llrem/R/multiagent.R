# Multi-agent REM-LLM simulation -----------------------------------------------

#' Compute the full hazard surface over all directed dyads
#'
#' Evaluates `exp(beta^T * x_ij)` for every directed dyad given the current
#' event history. Used internally by [run_multiagent_rem()] to derive sender
#' activation probabilities and the REM-hazard fallback for failed LLM parses.
#'
#' @param history     Data frame of past events (columns: time, sender,
#'   receiver, eveid).
#' @param base_riskset Precomputed base risk set (from [get_base_riskset()]).
#' @param beta        Named numeric vector of REM coefficients (from
#'   `coef(rem_fit)`).
#'
#' @return Data frame with columns `sender`, `receiver`, `hazard`.
#' @export
compute_hazard_surface <- function(history, base_riskset, beta,
                                    covariates = NULL) {
  if (nrow(history) == 0L) {
    return(data.frame(sender   = base_riskset$sender,
                      receiver = base_riskset$receiver,
                      hazard   = 1,
                      stringsAsFactors = FALSE))
  }
  dummy_ev <- data.frame(
    time = max(history$time) + 1L, sender = -1L, receiver = -1L,
    eveid = "-1_-1", stringsAsFactors = FALSE
  )
  rs <- riskset_cov2(history, dummy_ev, n_nonevents = NULL,
                     base_riskset = base_riskset, covariates = covariates)
  for (i in seq_len(ncol(rs)))
    if (is.numeric(rs[[i]])) { rs[[i]][is.na(rs[[i]])] <- 0; rs[[i]][!is.finite(rs[[i]])] <- 0 }
  beta_ok <- beta[!is.na(beta)]
  cn <- intersect(names(beta_ok), colnames(rs))
  X  <- as.matrix(rs[, cn, drop = FALSE])
  lp <- if (length(cn) > 0L) drop(X %*% beta_ok[cn]) else rep(0, nrow(rs))
  data.frame(sender   = rs$sender,
             receiver = rs$receiver,
             hazard   = exp(lp),
             stringsAsFactors = FALSE)
}

# query_agent ------------------------------------------------------------------
# Ask one LLM agent to choose a receiver.  On parse failure, falls back to
# sampling from the sender's REM hazard row rather than uniform random.

query_agent <- function(sender_id, valid_targets, history, cfg, hazard_row) {
  target_str <- paste(
    paste0(node_label(valid_targets, cfg), " (", valid_targets, ")"),
    collapse = ", "
  )
  msgs <- make_messages(
    paste0(cfg$context,
           "\n\nYou are roleplaying as one actor in this network. ",
           "Given the event history, decide who you will target next. ",
           "Reply with ONLY the integer ID of your target — no other text."),
    paste0("=== EVENT HISTORY ===\n", format_history(history, cfg),
           "\n\n=== YOUR TURN ===\n",
           "You are ", node_label(sender_id, cfg), " (ID: ", sender_id, ").\n",
           "Valid targets: ", target_str, "\n",
           "Which actor do you target? Reply with the integer ID only.")
  )
  raw  <- tryCatch(llm_call(cfg$provider, msgs, max_tokens = 15L),
                   error = function(e) NA_character_)
  recv <- first_valid_int(raw, valid_targets)
  if (is.na(recv)) recv <- label_match(raw, valid_targets, cfg)
  if (is.na(recv)) {
    h <- hazard_row
    h[is.na(h) | !is.finite(h)] <- 0
    if (sum(h) > 0) {
      recv <- sample(valid_targets, 1L, prob = h / sum(h))
    } else {
      recv <- sample(valid_targets, 1L)
    }
  }
  recv
}

#' Run a multi-agent REM-LLM conflict simulation
#'
#' At each tick, the REM hazard surface determines who acts (sender activation),
#' and individual LLM agents decide against whom (receiver selection). The two
#' components are complementary: the REM enforces endogenous structure; the LLM
#' contributes identity-conditioned content.
#'
#' @param edgelist        Prepared edgelist — used to infer the node set and
#'   optionally as the seed event.
#' @param rem_fit         Fitted `coxph` object (from [fit_rem()]).
#' @param cfg             A single [rem_cfg()] shared across all agents.
#'   Sender identity is injected dynamically at each query.
#' @param n_ticks         Number of simulation ticks.
#' @param lambda          Mean events per tick (Poisson rate). Estimated from
#'   `edgelist` inter-event times if `NULL`.
#' @param base_riskset    Precomputed base risk set. Built from `edgelist` if
#'   `NULL`.
#' @param cold_start      If `TRUE`, begin with a uniform hazard surface rather
#'   than seeding from the first empirical event. Cold-start runs have a burn-in
#'   period before endogenous structure accumulates.
#' @param seed            RNG seed.
#' @param checkpoint_path File path for incremental saves every 10 ticks.
#'   `NULL` disables checkpointing. Existing checkpoint is resumed automatically.
#'
#' @return A named list:
#'   \describe{
#'     \item{`sequence`}{Data frame of generated events (time, sender, receiver,
#'       eveid).}
#'     \item{`sender_probs`}{`n_ticks` × `N` matrix of activation probabilities
#'       \eqn{p_i(t)} at each tick.}
#'     \item{`n_events_per_tick`}{Integer vector of \eqn{K_t} values.}
#'   }
#' @export
run_multiagent_rem <- function(edgelist,
                                rem_fit,
                                cfg,
                                n_ticks,
                                lambda          = NULL,
                                base_riskset    = NULL,
                                cold_start      = FALSE,
                                seed            = 8214L,
                                checkpoint_path = NULL,
                                covariates      = NULL) {
  all_nodes    <- get_nodes(edgelist, NULL)
  base_riskset <- get_base_riskset(base_riskset, all_nodes)
  beta         <- coef(rem_fit)

  if (is.null(lambda)) {
    tr <- diff(range(edgelist$time))
    lambda <- if (tr > 0) nrow(edgelist) / tr else 1
  }

  set.seed(seed)

  # Initialise history and tick clock
  if (cold_start) {
    history   <- edgelist[0L, c("time", "sender", "receiver", "eveid")]
    tick_time <- 1L
  } else {
    history   <- edgelist[1L, c("time", "sender", "receiver", "eveid")]
    tick_time <- history$time[1L] + 1L
  }

  sequence_rows     <- list()
  sender_probs_list <- vector("list", n_ticks)
  n_events_per_tick <- integer(n_ticks)
  start_tick        <- 1L

  # Checkpoint recovery
  if (!is.null(checkpoint_path) && file.exists(checkpoint_path)) {
    env <- new.env(parent = emptyenv())
    load(checkpoint_path, envir = env)
    history           <- env$history
    sequence_rows     <- env$sequence_rows
    sender_probs_list <- env$sender_probs_list
    n_events_per_tick <- env$n_events_per_tick
    start_tick        <- env$next_tick
    tick_time         <- env$tick_time
    message("[llrem] Resuming multi-agent run from tick ", start_tick)
  }

  pb <- utils::txtProgressBar(min = start_tick, max = n_ticks, style = 3)

  for (tick in start_tick:n_ticks) {
    # Step 1 & 2: hazard surface → sender activation probabilities
    hz  <- compute_hazard_surface(history, base_riskset, beta, covariates)
    H_i <- tapply(hz$hazard, hz$sender, sum)
    H_i <- H_i[as.character(all_nodes)]
    H_i[is.na(H_i)] <- 0
    total_H <- sum(H_i)
    p_i <- if (total_H > 0) H_i / total_H else rep(1 / length(all_nodes), length(all_nodes))
    sender_probs_list[[tick]] <- as.numeric(p_i)

    # Step 3: number of events this tick
    K_t <- stats::rpois(1L, lambda)
    if (K_t == 0L) {
      n_events_per_tick[tick] <- 0L
      tick_time <- tick_time + 1L
      utils::setTxtProgressBar(pb, tick)
      next
    }
    K_t <- min(K_t, length(all_nodes))
    n_events_per_tick[tick] <- K_t

    # Step 4: sample K_t senders without replacement
    sampled_senders <- sample(all_nodes, size = K_t, replace = FALSE, prob = p_i)

    # Step 5: query each agent (all see the same pre-tick history)
    tick_rows <- vector("list", K_t)
    for (s_idx in seq_len(K_t)) {
      sender_id     <- sampled_senders[s_idx]
      valid_targets <- setdiff(all_nodes, sender_id)
      hz_sub        <- hz[hz$sender == sender_id, ]
      hz_row        <- hz_sub$hazard[match(valid_targets, hz_sub$receiver)]
      hz_row[is.na(hz_row)] <- 0

      recv <- query_agent(sender_id, valid_targets, history, cfg,
                          hazard_row = hz_row)
      tick_rows[[s_idx]] <- data.frame(
        time = tick_time, sender = sender_id, receiver = recv,
        eveid = paste0(sender_id, "_", recv), stringsAsFactors = FALSE
      )
    }

    # Step 6: append all tick events (within-tick events invisible to co-senders)
    tick_df       <- do.call(rbind, tick_rows)
    sequence_rows <- c(sequence_rows, tick_rows)
    history       <- rbind(history, tick_df)
    tick_time     <- tick_time + 1L

    if (!is.null(checkpoint_path) && (tick %% 10L == 0L || tick == n_ticks)) {
      next_tick <- tick + 1L
      save(history, sequence_rows, sender_probs_list, n_events_per_tick,
           next_tick, tick_time, file = checkpoint_path)
    }

    utils::setTxtProgressBar(pb, tick)
  }
  close(pb)

  sequence <- do.call(rbind, sequence_rows)
  if (!cold_start && nrow(edgelist) >= 1L) {
    sequence <- rbind(edgelist[1L, c("time", "sender", "receiver", "eveid")],
                      sequence)
  }
  row.names(sequence) <- NULL

  sender_probs <- do.call(rbind, lapply(sender_probs_list, function(p) {
    if (is.null(p)) rep(NA_real_, length(all_nodes)) else p
  }))
  colnames(sender_probs) <- as.character(all_nodes)

  list(sequence          = sequence,
       sender_probs      = sender_probs,
       n_events_per_tick = n_events_per_tick)
}
