# Prompt utilities (internal) -------------------------------------------------

node_label <- function(id, cfg) {
  if (!is.null(cfg$node_labels)) {
    lbl <- cfg$node_labels[as.character(id)]
    if (!is.na(lbl)) return(lbl)
  }
  paste(cfg$actor_noun, id)
}

format_history <- function(past_events, cfg) {
  if (nrow(past_events) == 0L) return("(No prior events.)")
  window <- cfg$history_window
  if (!is.null(window) && nrow(past_events) > window)
    past_events <- utils::tail(past_events, window)
  lines <- paste0(
    "  Event ", past_events$time, ": ",
    node_label(past_events$sender,   cfg), " (", past_events$sender,   ") ",
    cfg$action_verb, " ",
    node_label(past_events$receiver, cfg), " (", past_events$receiver, ")."
  )
  paste(lines, collapse = "\n")
}

make_messages <- function(system_msg, user_msg) {
  list(list(role = "system", content = system_msg),
       list(role = "user",   content = user_msg))
}

first_valid_int <- function(text, valid_set) {
  if (length(text) == 0L || is.na(text) || !nzchar(text)) return(NA_integer_)
  nums <- suppressWarnings(
    as.integer(regmatches(text, gregexpr("-?[0-9]+", text))[[1L]])
  )
  nums <- nums[nums %in% valid_set]
  if (length(nums) > 0L) nums[1L] else NA_integer_
}

label_match <- function(text, valid_set, cfg) {
  for (v in valid_set)
    if (grepl(node_label(v, cfg), text, ignore.case = TRUE)) return(v)
  NA_integer_
}

# Query factories -------------------------------------------------------------

#' Query factory: LLM roleplays as the empirical sender
#'
#' Returns a `query_fn` closure in which the LLM is told it IS the empirical
#' sender at time *t* and asked to choose a receiver. The generated sequence
#' preserves the empirical sender ordering but uses LLM-selected receivers.
#'
#' @param cfg A `rem_cfg` object.
#'
#' @return A function with signature
#'   `f(past_events, empirical_event, all_nodes) -> c(sender, receiver)`.
#' @export
make_query_roleplay_sender <- function(cfg) {
  force(cfg)
  function(past_events, empirical_event, all_nodes) {
    sender_id     <- empirical_event$sender
    valid_targets <- setdiff(all_nodes, sender_id)
    target_str    <- paste(
      paste0(node_label(valid_targets, cfg), " (", valid_targets, ")"),
      collapse = ", "
    )
    msgs <- make_messages(
      paste0(cfg$context,
             "\n\nYou are roleplaying as one actor in this network. ",
             "Given the event history, decide who you will target next. ",
             "Reply with ONLY the integer ID of your target — no other text."),
      paste0("=== EVENT HISTORY ===\n", format_history(past_events, cfg),
             "\n\n=== YOUR TURN ===\n",
             "You are ", node_label(sender_id, cfg), " (ID: ", sender_id, ").\n",
             "Valid targets: ", target_str, "\n",
             "Which actor do you target? Reply with the integer ID only.")
    )
    raw  <- llm_call(cfg$provider, msgs, max_tokens = 15L)
    recv <- first_valid_int(raw, valid_targets)
    if (is.na(recv)) recv <- label_match(raw, valid_targets, cfg)
    if (is.na(recv)) {
      warning("[llrem] roleplay_sender: could not parse '", raw,
              "'; substituting random target.", call. = FALSE)
      recv <- sample(valid_targets, 1L)
    }
    c(sender = sender_id, receiver = recv)
  }
}

#' Query factory: LLM predicts the next event
#'
#' Returns a `query_fn` closure that shows the LLM the event history and asks
#' it to predict both the sender and receiver of the next event. Tests the
#' model's implicit understanding of network dynamics rather than its own
#' behavioural tendencies.
#'
#' @param cfg A `rem_cfg` object.
#'
#' @return A function with signature
#'   `f(past_events, empirical_event, all_nodes) -> c(sender, receiver)`.
#' @export
make_query_predict <- function(cfg) {
  force(cfg)
  function(past_events, empirical_event, all_nodes) {
    node_str <- paste(
      paste0(node_label(all_nodes, cfg), " (", all_nodes, ")"), collapse = ", "
    )
    msgs <- make_messages(
      paste0(cfg$context,
             "\n\nObserve the sequence of past events and predict the next one. ",
             "Reply with EXACTLY two integers separated by a single space: ",
             "SENDER_ID RECEIVER_ID. No other text."),
      paste0("=== EVENT HISTORY ===\n", format_history(past_events, cfg),
             "\n\n=== PREDICTION TASK ===\n",
             "Actors: ", node_str, "\n",
             "Who acts next, and who do they target? Reply: SENDER_ID RECEIVER_ID")
    )
    raw  <- llm_call(cfg$provider, msgs, max_tokens = 25L)
    nums <- if (length(raw) == 0L || is.na(raw) || !nzchar(raw)) integer(0L) else suppressWarnings(
      as.integer(regmatches(raw, gregexpr("-?[0-9]+", raw))[[1L]])
    )
    nums     <- nums[nums %in% all_nodes]
    sender   <- if (length(nums) >= 1L) nums[1L] else label_match(raw, all_nodes, cfg)
    receiver <- if (length(nums) >= 2L) nums[2L] else NA_integer_
    valid_recv <- setdiff(all_nodes, sender)
    if (is.na(receiver) || identical(receiver, sender))
      receiver <- label_match(raw, valid_recv, cfg)
    if (is.na(sender) || is.na(receiver) || sender == receiver) {
      warning("[llrem] predict: could not parse '", raw,
              "'; substituting random event.", call. = FALSE)
      sender   <- sample(all_nodes, 1L)
      receiver <- sample(setdiff(all_nodes, sender), 1L)
    }
    c(sender = sender, receiver = receiver)
  }
}

#' Query factory: LLM roleplays as a randomly assigned node
#'
#' Returns a `query_fn` closure that assigns a uniformly random node identity
#' at each step and asks the LLM to choose a receiver. Provides a baseline
#' that isolates how much the assigned identity shapes LLM choices, independent
#' of the empirical sender sequence.
#'
#' @param cfg A `rem_cfg` object.
#'
#' @return A function with signature
#'   `f(past_events, empirical_event, all_nodes) -> c(sender, receiver)`.
#' @export
make_query_roleplay_random <- function(cfg) {
  force(cfg)
  inner <- make_query_roleplay_sender(cfg)
  function(past_events, empirical_event, all_nodes) {
    random_event        <- empirical_event
    random_event$sender <- sample(all_nodes, 1L)
    inner(past_events, random_event, all_nodes)
  }
}
