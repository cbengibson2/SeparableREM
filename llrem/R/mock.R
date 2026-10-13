#' Create a mock LLM provider for testing
#'
#' Returns an `llm_provider` that calls a local `strategy_fn` instead of
#' making any network requests. Useful for testing pipeline logic and custom
#' `query_fn` closures without incurring API costs.
#'
#' Three ready-made strategies are provided:
#' - [mock_strategy_min_id()]: always picks the smallest valid target ID.
#' - [mock_strategy_max_id()]: always picks the largest valid target ID.
#' - [mock_strategy_highest_indegree()]: picks the valid target with the
#'   highest cumulative in-degree in the event history.
#'
#' @param strategy_fn Function with signature
#'   `f(messages, max_tokens, temperature) -> character`. The return value
#'   must be parseable by the corresponding query factory (a stringified
#'   integer for roleplay scenarios; two space-separated integers for the
#'   predict scenario).
#' @param name Display name for `print()` output.
#'
#' @return An `llm_provider` object.
#' @export
llm_provider_mock <- function(strategy_fn, name = "Mock") {
  force(strategy_fn)
  new_provider(name, "mock",
    \(messages, max_tokens = 15L, temperature = 0.7)
      strategy_fn(messages, max_tokens, temperature)
  )
}

#' Mock strategy: always pick the smallest valid target
#'
#' @return A strategy function for use with [llm_provider_mock()].
#' @export
mock_strategy_min_id <- function(messages, max_tokens, temperature) {
  ids <- parse_valid_ids(messages)
  if (length(ids) == 0L) return("1")
  if (is_predict_msg(messages)) {
    s <- min(ids); r <- min(ids[ids != s])
    return(paste(s, r))
  }
  as.character(min(ids))
}

#' Mock strategy: always pick the largest valid target
#'
#' @return A strategy function for use with [llm_provider_mock()].
#' @export
mock_strategy_max_id <- function(messages, max_tokens, temperature) {
  ids <- parse_valid_ids(messages)
  if (length(ids) == 0L) return("1")
  if (is_predict_msg(messages)) {
    s <- max(ids); r <- max(ids[ids != s])
    return(paste(s, r))
  }
  as.character(max(ids))
}

#' Mock strategy: pick the valid target with the highest in-degree
#'
#' Parses the event history embedded in the prompt to tally cumulative
#' in-degrees, then returns the valid target node with the highest count.
#' Falls back to the smallest valid ID when no history is available yet.
#'
#' @return A strategy function for use with [llm_provider_mock()].
#' @export
mock_strategy_highest_indegree <- function(messages, max_tokens, temperature) {
  ids   <- parse_valid_ids(messages)
  if (length(ids) == 0L) return("1")
  indeg <- parse_indegree(messages)

  pick_best <- function(candidates) {
    if (length(indeg) > 0L) {
      ranked <- as.integer(names(indeg))
      best   <- ranked[ranked %in% candidates]
      if (length(best) > 0L) return(best[1L])
    }
    min(candidates)
  }

  if (is_predict_msg(messages)) {
    sender   <- pick_best(ids)
    receiver <- pick_best(ids[ids != sender])
    paste(sender, receiver)
  } else {
    as.character(pick_best(ids))
  }
}

# Message parsing helpers (internal) ------------------------------------------

user_content <- function(messages) {
  for (m in messages) if (identical(m$role, "user")) return(m$content)
  ""
}

parens_ints <- function(text) {
  m <- regmatches(text, gregexpr("(?<=\\()\\d+(?=\\))", text, perl = TRUE))[[1L]]
  suppressWarnings(as.integer(m))
}

is_predict_msg <- function(messages) {
  grepl("SENDER_ID RECEIVER_ID", user_content(messages), fixed = TRUE)
}

parse_valid_ids <- function(messages) {
  lines <- strsplit(user_content(messages), "\n", fixed = TRUE)[[1L]]
  tline <- grep("Valid targets:|Actors:", lines, value = TRUE)
  if (length(tline) == 0L) return(integer(0L))
  parens_ints(tline[1L])
}

parse_indegree <- function(messages) {
  lines   <- strsplit(user_content(messages), "\n", fixed = TRUE)[[1L]]
  evlines <- grep("Event \\d+:", lines, value = TRUE)
  if (length(evlines) == 0L) return(integer(0L))
  recvs <- vapply(evlines, function(ln) {
    ids <- parens_ints(ln)
    if (length(ids) >= 2L) ids[2L] else NA_integer_
  }, integer(1L), USE.NAMES = FALSE)
  sort(table(recvs[!is.na(recvs)]), decreasing = TRUE)
}
