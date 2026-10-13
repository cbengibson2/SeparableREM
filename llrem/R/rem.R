# Conditional logit (discrete-choice REM) -------------------------------------
#
# For each event t the risk set is {occurred=1 dyad} U {sampled non-events}.
# The log-likelihood is  Σ_t  β'x_{i_t j_t}  -  log Σ_{(i,j) ∈ R_t} exp(β'x_{ij})
# which is exactly McFadden's conditional logit with one chosen alternative per
# choice situation.  We maximise via optim(method="BFGS") with analytic gradient.

.rem_clogit_negll <- function(beta, X, strata) {
  eta    <- drop(X %*% beta)               # linear predictor, length = nrow(X)
  # subtract per-stratum max for numerical stability
  s_max  <- tapply(eta, strata, max)[strata]
  eta_c  <- eta - s_max
  exp_e  <- exp(eta_c)
  log_denom <- log(tapply(exp_e, strata, sum))[strata]
  lp     <- eta_c - log_denom              # log P(choose this row | stratum)
  ll     <- sum(lp[X[, 1L] == X[, 1L]])   # dummy; will use occurred mask below
  ll
}

# Fit conditional logit REM via BFGS.
# rs      : risk-set data frame with column 'occurred' and covariate columns
# formula : one-sided ~ covariate formula (no Surv(), no strata())
# Returns a list mimicking enough of coxph to be useful: coef, vcov, loglik
.fit_clogit <- function(rs, formula) {
  covars <- all.vars(formula)
  X      <- as.matrix(rs[, covars, drop = FALSE])
  y      <- rs$occurred
  g      <- rs$time2   # stratum = event time (unique per event)

  # negative log-likelihood
  negll <- function(beta) {
    eta   <- drop(X %*% beta)
    s_max <- tapply(eta, g, max)[as.character(g)]
    ec    <- exp(eta - s_max)
    ld    <- log(tapply(ec, g, sum))[as.character(g)]
    -sum((eta - s_max - ld)[y == 1L])
  }

  # analytic gradient
  grad <- function(beta) {
    eta   <- drop(X %*% beta)
    s_max <- tapply(eta, g, max)[as.character(g)]
    ec    <- exp(eta - s_max)
    denom <- tapply(ec, g, sum)[as.character(g)]
    p     <- ec / denom                   # P(row | stratum)
    # gradient = Σ_t [ x_{chosen_t} - E[x | stratum t] ]
    # = X'(y - p)  (where y is the 0/1 occurred indicator)
    -drop(t(X) %*% (y - p))
  }

  beta0 <- rep(0, ncol(X))
  fit   <- stats::optim(beta0, negll, grad, method = "BFGS",
                        control = list(maxit = 500L),
                        hessian = TRUE)

  coef <- fit$par
  names(coef) <- covars

  # vcov from observed Fisher information
  vcov <- tryCatch(solve(fit$hessian), error = function(e) {
    matrix(NA_real_, ncol(X), ncol(X),
           dimnames = list(covars, covars))
  })
  dimnames(vcov) <- list(covars, covars)

  structure(
    list(coef   = coef,
         vcov   = vcov,
         loglik = c(negll(rep(0, ncol(X))), -fit$value),
         converged = fit$convergence == 0L,
         method = "clogit"),
    class = "rem_clogit"
  )
}

#' @export
print.rem_clogit <- function(x, ...) {
  se  <- sqrt(diag(x$vcov))
  z   <- x$coef / se
  tab <- data.frame(coef = x$coef, `se(coef)` = se, z = z,
                    `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)),
                    check.names = FALSE)
  cat("REM conditional logit\n\n")
  print(round(tab, 4L))
  if (!x$converged) cat("\n[Warning: optimizer did not converge]\n")
  invisible(x)
}

# covariate formula from a Surv formula (strips Surv() and strata())
.surv_to_covar_formula <- function(formula) {
  tt  <- terms(formula, specials = c("strata", "survival::strata"))
  labs <- attr(tt, "term.labels")
  # drop any term whose label contains "strata("
  labs <- labs[!grepl("strata\\(", labs)]
  stats::as.formula(paste("~", paste(labs, collapse = " + ")))
}

#' Compute sufficient-statistic means for a sequence without fitting a REM
#'
#' Builds covariate values only for the realized events in `edgelist`
#' (n_nonevents = 0L) and returns their column means.  Much faster than
#' [fit_rem()] since no Cox model is fitted and the risk set contains exactly
#' one row per event.  Use to measure convergence in sufficient-statistic space
#' rather than coefficient space.
#'
#' @param edgelist    Prepared edgelist (output of [prepare_edgelist()]).
#' @param rem_formula Survival formula — used to identify covariate columns.
#' @param base_riskset Precomputed base risk set; constructed from `edgelist`
#'   if `NULL`.
#'
#' @return Named numeric vector: mean of each REM covariate across events.
#' @export
compute_seq_stat_means <- function(edgelist, rem_formula, base_riskset = NULL,
                                    covariates = NULL) {
  all_nodes    <- get_nodes(edgelist, NULL)
  base_riskset <- get_base_riskset(base_riskset, all_nodes)

  pe <- remstats_generate(edgelist,
                           n            = length(all_nodes),
                           n_nonevents  = 0L,
                           base_riskset = base_riskset,
                           covariates   = covariates)
  pe <- add_survival_cols(pe)

  ss <- make_suffstats(rem_formula, pe, drop_strata = TRUE)
  colMeans(ss, na.rm = TRUE)
}

#' Fit a Relational Event Model to an edgelist
#'
#' Builds risk-set statistics for each event in the sequence, adds survival
#' model timing columns (`time1`, `time2`), and fits a Cox proportional hazards
#' model. All covariate computation is handled internally.
#'
#' @param edgelist Data frame with columns `time`, `sender`, `receiver`,
#'   `eveid`. Use [prepare_edgelist()] to create this from a raw matrix.
#' @param formula A `survival::Surv` formula. Must reference `time1` and
#'   `time2` (via `strata()`). Available covariates: `Imitation`,
#'   `Repetition`, `RecInitiation`, `RRecInitiation`, `Recencylog`, `RSndSnd`,
#'   `WeightedRec`, `SenderRecencyLog`, `NODSnd`, `NODRec`, `NIDSnd`, `NIDRec`,
#'   `IOviolation`, `NTriadic`, `NOTPSnd`, `NITPSnd`, `NOSPSnd`, `NISPSnd`.
#' @param n Number of nodes. Inferred from `edgelist` if `NULL`.
#' @param n_nonevents Sampled non-event risk-set size per event. `NULL` uses
#'   the full risk set (slow for large networks; 200 is a sensible default).
#' @param base_riskset Precomputed base risk-set data frame. Constructed from
#'   `n` if `NULL`. Pass a cached object when calling `fit_rem()` and
#'   [run_llm_rem()] in the same session to avoid redundant work.
#' @param cl Optional PSOCK cluster from `parallel::makeCluster()` for
#'   parallelised risk-set computation.
#' @param method Estimation method: `"cox"` (default) fits a Cox proportional
#'   hazards model via [survival::coxph()] with one stratum per event.
#'   `"clogit"` fits a conditional logit (discrete-choice) model via BFGS,
#'   which is substantially faster when the number of events is large (>500)
#'   because it avoids the per-stratum partial-likelihood factorisation.
#'   Coefficients are on the same scale and asymptotically equivalent under
#'   the ordinal likelihood; `"clogit"` does not use event timing.
#'
#' @return For `method = "cox"`: a `coxph` object.  For `method = "clogit"`:
#'   a `rem_clogit` object with elements `coef`, `vcov`, `loglik`,
#'   and `converged`.
#' @export
fit_rem <- function(edgelist,
                    formula,
                    n            = NULL,
                    n_nonevents  = 200L,
                    base_riskset = NULL,
                    cl           = NULL,
                    method       = c("cox", "clogit"),
                    covariates   = NULL) {
  method       <- match.arg(method)
  all_nodes    <- get_nodes(edgelist, n)
  base_riskset <- get_base_riskset(base_riskset, all_nodes)
  pe           <- remstats_generate(edgelist,
                                     n            = length(all_nodes),
                                     n_nonevents  = n_nonevents,
                                     cl           = cl,
                                     base_riskset = base_riskset,
                                     covariates   = covariates)
  pe <- add_survival_cols(pe)

  if (method == "clogit") {
    covar_formula <- .surv_to_covar_formula(formula)
    return(.fit_clogit(pe, covar_formula))
  }
  survival::coxph(formula, data = pe)
}

#' Generate an LLM event sequence and fit a REM
#'
#' At each time step *t* (starting from 2), `query_fn` is called with the
#' growing LLM-generated history, the empirical event at *t*, and the full
#' node set. The returned sender/receiver pair is appended to the sequence.
#' The sequence is seeded with the first empirical event so the LLM starts
#' with at least one data point. After generation, risk-set statistics are
#' built and a REM is fitted.
#'
#' @param edgelist Prepared edgelist (output of [prepare_edgelist()]).
#' @param formula Survival formula passed to [fit_rem()].
#' @param query_fn A closure produced by [make_query_roleplay_sender()],
#'   [make_query_predict()], or [make_query_roleplay_random()], or any function
#'   with the signature:
#'   `f(past_events, empirical_event, all_nodes) -> c(sender, receiver)`.
#' @param n_events Number of events to generate. `NULL` uses the full edgelist.
#' @param n Number of nodes. Inferred if `NULL`.
#' @param n_nonevents Sampled risk-set size passed to [fit_rem()].
#' @param base_riskset Precomputed base risk set.
#' @param cl Optional PSOCK cluster for parallelised risk-set computation.
#' @param seed RNG seed for random fallbacks and random-node sampling.
#' @param checkpoint_path File path for incremental saves every 10 events.
#'   `NULL` disables checkpointing.
#' @param use_empirical_history If `FALSE` (default), the LLM is shown its own
#'   previously generated events as history, so any divergence from reality
#'   compounds over time. If `TRUE`, the LLM is shown the empirical event
#'   history at each step, isolating the effect of the generation strategy from
#'   the effect of accumulated sequence drift.
#' @param rate_limit_secs Seconds to pause between LLM calls.
#' @param fit If `TRUE` (default), fits a REM on the generated sequence before
#'   returning. Set to `FALSE` to skip model fitting and return only the
#'   sequence and raw responses — useful when running many sequences and
#'   deferring all model fitting to a separate step.
#' @param method Passed to [fit_rem()]: `"cox"` (default) or `"clogit"`.
#'
#' @return A named list with elements:
#'   - `sequence`: generated event sequence as a data frame.
#'   - `responses`: raw LLM response text and parsed sender/receiver per event.
#'   - `rem`: fitted model (only present when `fit = TRUE`).
#' @export
run_llm_rem <- function(edgelist,
                         formula,
                         query_fn,
                         n_events              = NULL,
                         n                     = NULL,
                         n_nonevents           = 200L,
                         base_riskset          = NULL,
                         cl                    = NULL,
                         seed                  = 8214L,
                         checkpoint_path       = NULL,
                         use_empirical_history = FALSE,
                         rate_limit_secs       = 0.3,
                         fit                   = TRUE,
                         method                = c("cox", "clogit"),
                         covariates            = NULL) {
  method <- match.arg(method)
  all_nodes    <- get_nodes(edgelist, n)
  base_riskset <- get_base_riskset(base_riskset, all_nodes)

  n_ev <- if (is.null(n_events)) nrow(edgelist) else min(n_events, nrow(edgelist))
  emp  <- edgelist[seq_len(n_ev), ]

  set.seed(seed)
  llm_rows <- list()
  start_t  <- 2L

  if (!is.null(checkpoint_path) && file.exists(checkpoint_path)) {
    env <- new.env(parent = emptyenv())
    load(checkpoint_path, envir = env)
    llm_rows <- env$llm_rows
    start_t  <- length(llm_rows) + 2L
    message("[llrem] Resuming from checkpoint at event ", start_t)
  }

  llm_seq <- emp[1L, c("time", "sender", "receiver", "eveid")]

  if (start_t > n_ev) {
    message("[llrem] Checkpoint is complete; skipping generation.")
  } else {
    pb <- utils::txtProgressBar(min = start_t, max = n_ev, style = 3)

    for (t in start_t:n_ev) {
      history <- if (use_empirical_history) emp[seq_len(t - 1L), ] else llm_seq
      result <- query_fn(
        past_events     = history,
        empirical_event = emp[t, ],
        all_nodes       = all_nodes
      )
      new_row <- data.frame(
        time     = emp$time[t],
        sender   = as.integer(result[["sender"]]),
        receiver = as.integer(result[["receiver"]]),
        eveid    = paste0(result[["sender"]], "_", result[["receiver"]]),
        stringsAsFactors = FALSE
      )
      llm_rows[[t - 1L]] <- new_row
      llm_seq <- rbind(llm_seq, new_row[, c("time", "sender", "receiver", "eveid")])

      if (!is.null(checkpoint_path) && (t %% 10L == 0L || t == n_ev))
        save(llm_rows, file = checkpoint_path)

      utils::setTxtProgressBar(pb, t)
      Sys.sleep(rate_limit_secs)
    }
    close(pb)
  }

  sequence <- rbind(emp[1L, c("time", "sender", "receiver", "eveid")],
                    do.call(rbind, llm_rows))
  row.names(sequence) <- NULL

  out <- list(sequence  = sequence,
              responses = do.call(rbind, llm_rows))

  if (fit) {
    message("\n[llrem] Fitting REM on generated sequence...")
    out$rem <- fit_rem(sequence, formula,
                       n = length(all_nodes), n_nonevents = n_nonevents,
                       base_riskset = base_riskset, cl = cl, method = method,
                       covariates = covariates)
  }

  out
}
