#' Extract sufficient-statistic means from a REM risk-set data frame
#'
#' Returns a data frame with one row per value of `occurred` (0/1), containing
#' the mean of each REM covariate.  Used as input to [make_behavior_queries()].
#'
#' @param rem_formula Survival formula whose RHS terms are the covariates.
#' @param data        Risk-set data frame (expanded counting-process format)
#'   with an `occurred` column.
#' @param drop_strata If `TRUE` (default), drop strata columns from the model
#'   matrix before computing means.
#'
#' @return A data frame with columns `occurred` plus one column per covariate.
#' @export
make_suffstats <- function(rem_formula, data, drop_strata = TRUE) {
  tt <- terms(rem_formula, specials = "strata", data = data)
  X  <- model.matrix(tt, data = data)

  if ("(Intercept)" %in% colnames(X))
    X <- X[, colnames(X) != "(Intercept)", drop = FALSE]

  if (drop_strata)
    X <- X[, !grepl("^strata\\(", colnames(X)), drop = FALSE]

  as.data.frame(X)
}

#' Compute per-term mean sufficient statistics for realized events
#'
#' Aggregates a risk-set data frame by `occurred` and returns means of the REM
#' covariates, plus the median of `Recency` for realized events (if present).
#'
#' @param rem_formula Survival formula.
#' @param poss_edges  Risk-set data frame containing `occurred` and `Recency`.
#'
#' @return A data frame suitable for passing to [make_behavior_queries()].
#' @export
make_stat_means <- function(rem_formula, poss_edges) {
  ss <- make_suffstats(rem_formula, poss_edges)
  sm <- aggregate(
    . ~ occurred,
    data = cbind(occurred = poss_edges$occurred, ss),
    FUN  = mean
  )
  if ("Recency" %in% names(poss_edges)) {
    sm$RecencyMedian <- c(
      median(poss_edges$Recency[poss_edges$occurred == 0]),
      median(poss_edges$Recency[poss_edges$occurred == 1 &
                                  poss_edges$RecInitiation == 1])
    )
  }
  sm
}

#' Generate plain-language behavioral guidelines from empirical sufficient stats
#'
#' Translates the mean sufficient statistics for realized conflict events into
#' numbered natural-language instructions suitable for inclusion in an LLM
#' system prompt.  Intended to ground LLM predictions in the empirical base
#' rates of each relational pattern.
#'
#' @param stat_means  Data frame returned by [make_stat_means()]: one row per
#'   value of `occurred`, one column per REM term.
#' @param event_value The value of `occurred` that identifies realized events
#'   (default 1).
#' @param digits      Decimal digits for percentage rounding (default 0).
#'
#' @return A single character string: a numbered list of behavioral guidelines.
#' @export
make_behavior_queries <- function(stat_means, event_value = 1, digits = 0) {
  event_row          <- stat_means[stat_means$occurred == event_value, ,
                                   drop = FALSE]
  event_row$occurred <- NULL

  pct <- function(x) paste0(round(100 * x, digits), "%")
  out <- character(0L)

  if ("Repetition" %in% names(event_row))
    out <- c(out, sprintf(
      "About %s of realized events from i to j had at least one prior i \u2192 j event. Prefer previously used i \u2192 j dyads over entirely new dyads.",
      pct(event_row$Repetition)
    ))

  if ("RecInitiation" %in% names(event_row))
    out <- c(out, sprintf(
      "About %s of realized events from i to j had at least one prior j \u2192 i event. Prefer targets that have previously attacked i.",
      pct(event_row$RecInitiation)
    ))

  if ("IOviolation" %in% names(event_row))
    out <- c(out, sprintf(
      "Only about %s of realized events are interaction-order violations — events where a new actor targets the previous sender (AB-XA), the sender redirects to a new target (AB-XY),
      or a new actor targets the previous receiver (PSAB-XB). In this data, these participation shift patterns are strongly suppressed: after A attacks B, the sequence ten
ds NOT to jump to a new actor or a new target pairing. Reduce how often you introduce new participants or redirect aggression. For example: 69 attacked 10, but shortly after 38 attack
ed 69 (a new actor turning on the previous sender — PSAB-XA). These participation shifts are less common than simple repetition.",
      pct(event_row$IOviolation)
    ))

  if ("Imitation" %in% names(event_row))
    out <- c(out, sprintf(
      "About %s of realized events involve imitation. Do not often choose a target simply because others recently attacked that target.",
      pct(event_row$Imitation)
    ))

  if ("NTriadic" %in% names(event_row))
    out <- c(out, sprintf(
      "About %s of realized events involve the triadic pattern. Use triadic closure only rarely.",
      pct(event_row$NTriadic)
    ))

  if ("Recencylog" %in% names(event_row) && "RecencyMedian" %in% names(event_row))
    out <- c(out, sprintf(
      "Actors typically reciprocate an act within %s timepoints or less. Prefer i \u2192 j events that have had a recent j \u2192 i event.",
      event_row$RecencyMedian
    ))

  paste(seq_along(out), out, sep = ". ", collapse = "\n")
}
