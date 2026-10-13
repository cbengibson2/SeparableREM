# Exogenous covariate support for REM risk sets --------------------------------
#
# Four covariate types, each an S3 list with class c("<type>", "rem_cov"):
#   cov_node_static    – one value per node
#   cov_node_temporal  – one value per (node, time-window)
#   cov_dyad_static    – one value per directed dyad
#   cov_dyad_temporal  – one value per (dyad, time-window)
#
# All four are collected into a rem_covariates() container that validates
# names and passes as a single argument through the pipeline.

.reserved_cols <- c(
  "sender", "receiver", "eid", "occurred", "time", "time1", "time2", "eveid",
  "senOutDg", "dyad_occurrence", "Repetition", "Recency", "Recencylog",
  "PSABAB", "PSABBA", "PSABAY", "PSABXA", "PSABXB", "PSABBY", "PSABXY",
  "IOviolation", "recInDg", "PreImitation", "Imitation",
  "OSP", "ISP", "ITP", "OTP", "NOTPSnd", "NITPSnd", "NOSPSnd", "NISPSnd",
  "triadic_closure", "NTriadic", "RecInitiation", "RRecInitiation", "RSndSnd",
  "w_ji", "WeightedRec", "SenderRecencyLog",
  "NODSnd", "NODRec", "NIDSnd", "NIDRec", "senInDg", "recOutDg", "n_norm"
)

# Constructors -----------------------------------------------------------------

#' Create a static node-level covariate
#'
#' Adds one or two columns to each risk-set row based on the sender's and/or
#' receiver's value in a node attribute table.
#'
#' @param name Character scalar. Column name(s) added to the risk set.
#'   When `role = "both"`, two columns are added: `<name>_snd` and
#'   `<name>_rec`. Otherwise a single column named `<name>` is added.
#' @param data Data frame with columns `node` (integer-coercible) and
#'   `value` (numeric). One row per node.
#' @param role One of `"both"` (default), `"sender"`, or `"receiver"`.
#' @param warn_missing Logical. Warn when a node in the risk set has no
#'   matching row in `data` (default `TRUE`).
#'
#' @return An object of class `c("cov_node_static", "rem_cov")`.
#' @export
cov_node_static <- function(name, data, role = c("both", "sender", "receiver"),
                             warn_missing = TRUE) {
  role <- match.arg(role)
  stopifnot(is.character(name), length(name) == 1L)
  stopifnot(is.data.frame(data))
  if (!all(c("node", "value") %in% names(data)))
    stop("'data' must have columns 'node' and 'value'.", call. = FALSE)
  data$node <- as.integer(data$node)
  structure(
    list(name = name, data = data, role = role, warn_missing = warn_missing),
    class = c("cov_node_static", "rem_cov")
  )
}

#' Create a temporal node-level covariate
#'
#' Like [cov_node_static()] but values change at known time breakpoints.
#' The active value at event time `t` is the one with the largest
#' `time_start <= t` (step function).
#'
#' @param name Character scalar.
#' @param data Data frame with columns `node` (integer-coercible),
#'   `time_start` (numeric), and `value` (numeric). Rows need not be
#'   sorted; sorting is done internally.
#' @param role One of `"both"` (default), `"sender"`, or `"receiver"`.
#' @param warn_missing Logical.
#'
#' @return An object of class `c("cov_node_temporal", "rem_cov")`.
#' @export
cov_node_temporal <- function(name, data,
                               role = c("both", "sender", "receiver"),
                               warn_missing = TRUE) {
  role <- match.arg(role)
  stopifnot(is.character(name), length(name) == 1L)
  stopifnot(is.data.frame(data))
  if (!all(c("node", "time_start", "value") %in% names(data)))
    stop("'data' must have columns 'node', 'time_start', and 'value'.",
         call. = FALSE)
  data$node <- as.integer(data$node)
  # Pre-split by node for O(log W) lookup at merge time instead of O(m)
  data       <- data[order(data$node, data$time_start), ]
  by_node    <- split(data[, c("time_start", "value")], data$node)
  structure(
    list(name = name, data = data, by_node = by_node,
         role = role, warn_missing = warn_missing),
    class = c("cov_node_temporal", "rem_cov")
  )
}

#' Create a static dyad-level covariate
#'
#' Adds one column to each risk-set row based on the (sender, receiver) pair.
#'
#' @param name Character scalar. Column name added to the risk set.
#' @param data Data frame with columns `sender` (integer-coercible),
#'   `receiver` (integer-coercible), and `value` (numeric). Need not be a
#'   complete cross-product; unmatched dyads receive 0 (with optional warning).
#' @param warn_missing Logical.
#'
#' @return An object of class `c("cov_dyad_static", "rem_cov")`.
#' @export
cov_dyad_static <- function(name, data, warn_missing = TRUE) {
  stopifnot(is.character(name), length(name) == 1L)
  stopifnot(is.data.frame(data))
  if (!all(c("sender", "receiver", "value") %in% names(data)))
    stop("'data' must have columns 'sender', 'receiver', and 'value'.",
         call. = FALSE)
  data$sender   <- as.integer(data$sender)
  data$receiver <- as.integer(data$receiver)
  # Pre-compute eid key once at construction
  data$.eid <- paste0(data$sender, "_", data$receiver)
  structure(
    list(name = name, data = data, warn_missing = warn_missing),
    class = c("cov_dyad_static", "rem_cov")
  )
}

#' Create a temporal dyad-level covariate
#'
#' Like [cov_dyad_static()] but values change at known time breakpoints.
#'
#' @param name Character scalar.
#' @param data Data frame with columns `sender`, `receiver`, `time_start`,
#'   and `value`. Rows need not be sorted.
#' @param warn_missing Logical.
#'
#' @return An object of class `c("cov_dyad_temporal", "rem_cov")`.
#' @export
cov_dyad_temporal <- function(name, data, warn_missing = TRUE) {
  stopifnot(is.character(name), length(name) == 1L)
  stopifnot(is.data.frame(data))
  if (!all(c("sender", "receiver", "time_start", "value") %in% names(data)))
    stop("'data' must have columns 'sender', 'receiver', 'time_start', and 'value'.",
         call. = FALSE)
  data$sender   <- as.integer(data$sender)
  data$receiver <- as.integer(data$receiver)
  data$.eid     <- paste0(data$sender, "_", data$receiver)
  data          <- data[order(data$.eid, data$time_start), ]
  by_dyad       <- split(data[, c("time_start", "value")], data$.eid)
  structure(
    list(name = name, data = data, by_dyad = by_dyad, warn_missing = warn_missing),
    class = c("cov_dyad_temporal", "rem_cov")
  )
}

#' Collect covariate objects for use in fit_rem
#'
#' Validates that all covariate names are unique, do not collide with
#' endogenous REM column names, and are proper `rem_cov` objects.
#'
#' @param ... Any number of objects created by [cov_node_static()],
#'   [cov_node_temporal()], [cov_dyad_static()], or [cov_dyad_temporal()].
#'
#' @return An object of class `"rem_covariates"`.
#' @export
rem_covariates <- function(...) {
  covs <- list(...)
  if (length(covs) == 0L)
    return(structure(list(covariates = list()), class = "rem_covariates"))

  for (i in seq_along(covs)) {
    if (!inherits(covs[[i]], "rem_cov"))
      stop("All arguments to rem_covariates() must be rem_cov objects ",
           "(created by cov_node_static, cov_node_temporal, ",
           "cov_dyad_static, or cov_dyad_temporal).", call. = FALSE)
  }

  all_names <- vapply(covs, `[[`, character(1L), "name")

  # node covariates with role="both" produce two columns: name_snd, name_rec
  col_names <- unlist(lapply(covs, function(cv) {
    if (inherits(cv, c("cov_node_static", "cov_node_temporal")) &&
        cv$role == "both") {
      c(paste0(cv$name, "_snd"), paste0(cv$name, "_rec"))
    } else {
      cv$name
    }
  }))

  dupes <- col_names[duplicated(col_names)]
  if (length(dupes))
    stop("Duplicate covariate column name(s): ",
         paste(unique(dupes), collapse = ", "), call. = FALSE)

  all_to_check  <- unique(c(all_names, col_names))
  reserved_hits <- all_to_check[all_to_check %in% .reserved_cols]
  if (length(reserved_hits))
    stop("Covariate name(s) collide with reserved REM columns: ",
         paste(reserved_hits, collapse = ", "), call. = FALSE)

  structure(list(covariates = covs), class = "rem_covariates")
}

#' @export
print.rem_covariates <- function(x, ...) {
  covs <- x$covariates
  if (length(covs) == 0L) {
    cat("rem_covariates: (empty)\n")
    return(invisible(x))
  }
  cat("rem_covariates with", length(covs), "covariate(s):\n")
  for (cv in covs) {
    type <- sub("cov_", "", class(cv)[1L])
    extra <- if (inherits(cv, c("cov_node_static", "cov_node_temporal")))
      paste0(", role=", cv$role) else ""
    cat(sprintf("  %-20s  [%s%s]\n", cv$name, type, extra))
  }
  invisible(x)
}

# Internal merge worker --------------------------------------------------------

# Step-function lookup: index of last time_start <= t, or 0 if none.
.find_window <- function(t, time_starts) {
  findInterval(t, time_starts, rightmost.closed = FALSE, left.open = FALSE)
}

# Look up node values for a vector of node IDs given one covariate and event time.
.node_vals <- function(nodes, cv, t) {
  if (inherits(cv, "cov_node_static")) {
    idx <- match(nodes, cv$data$node)
    cv$data$value[idx]
  } else {
    # temporal: per-node step-function lookup
    vals <- vapply(nodes, function(nd) {
      d <- cv$by_node[[as.character(nd)]]
      if (is.null(d)) return(NA_real_)
      w <- .find_window(t, d$time_start)
      if (w == 0L) NA_real_ else d$value[w]
    }, numeric(1L))
    vals
  }
}

# Append all covariates in a rem_covariates object to a riskset data frame.
# event_time: scalar integer/numeric clock time of the current event.
attach_covariates <- function(riskset, covariates, event_time) {
  if (is.null(covariates) || length(covariates$covariates) == 0L)
    return(riskset)

  for (cv in covariates$covariates) {

    if (inherits(cv, c("cov_node_static", "cov_node_temporal"))) {
      if (cv$role %in% c("sender", "both")) {
        vals <- .node_vals(riskset$sender, cv, event_time)
        n_miss <- sum(is.na(vals))
        if (n_miss > 0L && cv$warn_missing)
          warning(sprintf("cov '%s' (_snd): %d risk-set row(s) have no node match",
                          cv$name, n_miss), call. = FALSE)
        col <- if (cv$role == "both") paste0(cv$name, "_snd") else cv$name
        riskset[[col]] <- vals
      }
      if (cv$role %in% c("receiver", "both")) {
        vals <- .node_vals(riskset$receiver, cv, event_time)
        n_miss <- sum(is.na(vals))
        if (n_miss > 0L && cv$warn_missing)
          warning(sprintf("cov '%s' (_rec): %d risk-set row(s) have no node match",
                          cv$name, n_miss), call. = FALSE)
        col <- if (cv$role == "both") paste0(cv$name, "_rec") else cv$name
        riskset[[col]] <- vals
      }

    } else if (inherits(cv, "cov_dyad_static")) {
      idx  <- match(riskset$eid, cv$data$.eid)
      vals <- cv$data$value[idx]
      n_miss <- sum(is.na(vals))
      if (n_miss > 0L && cv$warn_missing)
        warning(sprintf("cov '%s': %d risk-set row(s) have no dyad match",
                        cv$name, n_miss), call. = FALSE)
      riskset[[cv$name]] <- vals

    } else if (inherits(cv, "cov_dyad_temporal")) {
      vals <- vapply(riskset$eid, function(eid) {
        d <- cv$by_dyad[[eid]]
        if (is.null(d)) return(NA_real_)
        w <- .find_window(event_time, d$time_start)
        if (w == 0L) NA_real_ else d$value[w]
      }, numeric(1L))
      n_miss <- sum(is.na(vals))
      if (n_miss > 0L && cv$warn_missing)
        warning(sprintf("cov '%s': %d risk-set row(s) have no dyad/time match",
                        cv$name, n_miss), call. = FALSE)
      riskset[[cv$name]] <- vals
    }
  }
  riskset
}

#' Attach exogenous covariates to a pre-built risk-set data frame
#'
#' Convenience wrapper for applying a [rem_covariates()] object to the output
#' of `remstats_generate()` after the fact.  Useful for testing or for
#' workflows that build the risk set once and experiment with different
#' covariate specifications.
#'
#' @param pe          Data frame returned by `remstats_generate()` (or
#'   equivalently, the `pe` object inside [fit_rem()]), with a `time` column.
#' @param covariates  A [rem_covariates()] object.
#'
#' @return `pe` with additional covariate columns appended.
#' @export
attach_covariates_df <- function(pe, covariates) {
  if (is.null(covariates) || length(covariates$covariates) == 0L)
    return(pe)
  chunks <- split(pe, pe$time)
  chunks <- lapply(chunks, function(chunk) {
    t <- chunk$time[1L]
    attach_covariates(chunk, covariates, t)
  })
  do.call(rbind, chunks)
}
