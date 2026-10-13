#' Normalise a raw edgelist to the standard three-column format
#'
#' Accepts a matrix or data frame with columns in order: time, sender,
#' receiver. Returns a data frame with named integer columns and an `eveid`
#' dyad key (`"sender_receiver"`).
#'
#' @param raw Matrix or data frame with at least 3 columns.
#'
#' @return A data frame with columns `time` (integer), `sender` (integer),
#'   `receiver` (integer), and `eveid` (character).
#' @export
prepare_edgelist <- function(raw) {
  df <- as.data.frame(raw)
  if (ncol(df) < 3L)
    stop("`raw` must have at least 3 columns: time, sender, receiver.",
         call. = FALSE)
  colnames(df)[1:3] <- c("time", "sender", "receiver")
  df$time     <- as.integer(df$time)
  df$sender   <- as.integer(df$sender)
  df$receiver <- as.integer(df$receiver)
  df$eveid    <- paste0(df$sender, "_", df$receiver)
  df
}

# Internal helpers ------------------------------------------------------------

get_nodes <- function(edgelist, n_hint = NULL) {
  nodes <- sort(unique(c(edgelist$sender, edgelist$receiver)))
  if (!is.null(n_hint) && n_hint > length(nodes)) nodes <- seq_len(n_hint)
  nodes
}

get_base_riskset <- function(base_riskset, all_nodes) {
  if (!is.null(base_riskset)) return(base_riskset)
  br     <- expand.grid(sender = all_nodes, receiver = all_nodes)
  br     <- br[br$sender != br$receiver, ]
  br$eid <- paste0(br$sender, "_", br$receiver)
  br
}

add_survival_cols <- function(pe) {
  pe$time1 <- 1L
  pe$time2 <- pe$time
  pe
}

# riskset_cov2 ----------------------------------------------------------------
# Compute REM covariates for one event given all past events.
#
# Covariates:
#   Repetition    dyad_occurrence / sender out-degree
#   Recency       time of last reverse-dyad event
#   Recencylog    -log1p(Recency)
#   IOviolation   interaction-order violation (Simmelian sequence)
#   Imitation     normalised indegree excluding own prior ties
#   OTP / ITP / OSP / ISP   raw two-path / shared-partner counts
#   NOTPSnd       OTP normalised by per-sender max (transitivity)
#   NITPSnd       ITP normalised by per-sender max (cyclic closure)
#   NOSPSnd       OSP normalised by per-sender max (outbound shared partners)
#   NISPSnd       ISP normalised by per-sender max (inbound shared partners)
#   NTriadic      normalised sum of all four (backward-compat composite)
#   NODSnd        normalised out-degree of sender:   senOutDg / (n-1)
#   NODRec        normalised out-degree of receiver: recOutDg / (n-1)
#   NIDSnd        normalised in-degree  of sender:   senInDg  / (n-1)
#   NIDRec        normalised in-degree  of receiver: recInDg  / (n-1)
#   RecInitiation has this dyad occurred before? (binary)
#   RRecInitiation has the reverse dyad ever occurred? (reciprocity binary)
#   RSndSnd       -log1p(time of last same-dyad event); 0 if never fired
#   WeightedRec      count of reverse-dyad events / per-sender max
#   SenderRecencyLog log1p(time_of_last_i→j * RecInitiation); 0 if never fired
#   occurred      1 for the focal event, 0 for non-events

riskset_cov2 <- function(past_events, current_event,
                          n_nonevents = NULL, base_riskset = NULL,
                          covariates = NULL) {
  riskset <- base_riskset
  ev      <- match(current_event$eveid, riskset$eid)

  if (!is.null(n_nonevents) && n_nonevents < nrow(riskset)) {
    keep    <- unique(c(ev, sample.int(nrow(riskset), n_nonevents)))
    riskset <- riskset[keep, ]
  }

  # Repetition
  outd                    <- table(past_events$sender)
  riskset$senOutDg        <- outd[match(riskset$sender, names(outd))]
  docc                    <- table(past_events$eveid)
  riskset$dyad_occurrence <- docc[match(riskset$eid, names(docc))]
  riskset$Repetition      <- riskset$dyad_occurrence / riskset$senOutDg

  # Recency
  last_time          <- tapply(past_events$time, past_events$eveid, max)
  recip_ids          <- paste0(riskset$receiver, "_", riskset$sender)
  riskset$Recency    <- last_time[match(recip_ids, names(last_time))]
  riskset$Recency[is.na(riskset$Recency)] <- 0
  riskset$Recencylog <- -log1p(riskset$Recency)

  # Interaction-order violation
  lastsender   <- utils::tail(past_events$sender,   1L)
  lastreceiver <- utils::tail(past_events$receiver, 1L)
  riskset$PSABAB <- riskset$sender   == lastsender   & riskset$receiver == lastreceiver
  riskset$PSABBA <- riskset$receiver == lastsender   & riskset$sender   == lastreceiver
  riskset$PSABAY <- riskset$sender   == lastsender   & !riskset$PSABAB
  riskset$PSABXA <- riskset$receiver == lastsender   & !riskset$PSABBA
  riskset$PSABXB <- riskset$receiver == lastreceiver & !riskset$PSABAB
  riskset$PSABBY <- riskset$sender   == lastreceiver & !riskset$PSABBA
  riskset$PSABXY <- (riskset$PSABXB + riskset$PSABXA + riskset$PSABAY +
                       riskset$PSABAB + riskset$PSABBA + riskset$PSABBY) == 0
  riskset$IOviolation <- riskset$PSABXA | riskset$PSABXY | riskset$PSABXB

  # Imitation
  indg                 <- table(past_events$receiver)
  riskset$recInDg      <- indg[match(riskset$receiver, names(indg))]
  riskset$PreImitation <- riskset$recInDg - riskset$dyad_occurrence
  riskset$PreImitation[is.na(riskset$PreImitation)] <- 0
  imimx                <- tapply(riskset$PreImitation, riskset$sender, max)
  imimx                <- imimx[match(riskset$sender, names(imimx))]
  riskset$Imitation    <- riskset$PreImitation / imimx

  # Triadic closure (OSP, ISP, ITP, OTP) via C++
  adj              <- build_adjacency(past_events$sender, past_events$receiver)
  counts           <- triadic_counts(adj, riskset$sender, riskset$receiver)
  riskset$OSP <- counts[, "OSP"]
  riskset$ISP <- counts[, "ISP"]
  riskset$ITP <- counts[, "ITP"]
  riskset$OTP <- counts[, "OTP"]

  # Per-sender-max normalisation (mirrors relevent's *Snd convention)
  norm_by_sender_max <- function(x, sender) {
    mx <- tapply(x, sender, max)
    mx <- mx[match(sender, names(mx))]
    ifelse(mx > 0, x / mx, 0)
  }
  riskset$NOTPSnd <- norm_by_sender_max(riskset$OTP, riskset$sender)
  riskset$NITPSnd <- norm_by_sender_max(riskset$ITP, riskset$sender)
  riskset$NOSPSnd <- norm_by_sender_max(riskset$OSP, riskset$sender)
  riskset$NISPSnd <- norm_by_sender_max(riskset$ISP, riskset$sender)

  riskset$triadic_closure <- riskset$OSP + riskset$ISP + riskset$ITP + riskset$OTP
  riskset$triadic_closure[is.na(riskset$triadic_closure)] <- 0L
  mx_tri           <- max(riskset$triadic_closure)
  riskset$NTriadic <- if (mx_tri > 0) riskset$triadic_closure / mx_tri else 0

  # Normalised degree statistics (match relevent NODSnd/NODRec/NIDSnd/NIDRec)
  # senOutDg and recInDg are already computed above; add the other two here.
  # n_norm = (n-1): derive from base_riskset node set, not just past_events,
  # so early events don't inflate the ratio above 1.
  n_norm                  <- max(length(unique(c(base_riskset$sender,
                                                 base_riskset$receiver))) - 1L,
                                 1L)
  senInDg                 <- table(past_events$receiver)
  recOutDg                <- table(past_events$sender)
  riskset$NODSnd <- as.numeric(riskset$senOutDg) / n_norm
  riskset$NODSnd[is.na(riskset$NODSnd)] <- 0
  riskset$NODRec <- as.numeric(recOutDg[match(riskset$receiver, names(recOutDg))]) / n_norm
  riskset$NODRec[is.na(riskset$NODRec)] <- 0
  riskset$NIDSnd <- as.numeric(senInDg[match(riskset$sender, names(senInDg))]) / n_norm
  riskset$NIDSnd[is.na(riskset$NIDSnd)] <- 0
  riskset$NIDRec <- as.numeric(riskset$recInDg) / n_norm
  riskset$NIDRec[is.na(riskset$NIDRec)] <- 0

  riskset$occurred        <- as.numeric(riskset$eid %in% current_event$eveid)
  riskset$time            <- current_event$time
  riskset$RecInitiation   <- as.numeric(riskset$dyad_occurrence > 0)
  riskset$RRecInitiation  <- as.numeric(!is.na(last_time[match(recip_ids, names(last_time))]))
  snd_last                <- last_time[match(riskset$eid, names(last_time))]
  snd_last[is.na(snd_last)] <- 0
  riskset$RSndSnd         <- -log1p(snd_last)

  # WeightedRec: count of reverse-dyad (j→i) events prior to t, normalised
  # by the per-sender maximum across the current risk set
  rev_docc                <- docc[match(recip_ids, names(docc))]
  rev_docc[is.na(rev_docc)] <- 0L
  riskset$w_ji            <- as.integer(rev_docc)
  riskset$WeightedRec      <- norm_by_sender_max(riskset$w_ji, riskset$sender)
  riskset$SenderRecencyLog <- log1p(snd_last * riskset$RecInitiation)

  if (!is.null(covariates))
    riskset <- attach_covariates(riskset, covariates, current_event$time)

  riskset
}

# remstats_generate -----------------------------------------------------------

remstats_generate <- function(edgelist, n = NULL, n_nonevents = NULL,
                               cl = NULL, base_riskset = NULL,
                               covariates = NULL) {
  if (is.null(n))
    n <- length(unique(c(edgelist$sender, edgelist$receiver)))
  if (is.null(n_nonevents))
    n_nonevents <- n * (n - 1L)

  if (is.null(base_riskset)) {
    base_riskset     <- expand.grid(sender = seq_len(n), receiver = seq_len(n))
    base_riskset     <- base_riskset[base_riskset$sender != base_riskset$receiver, ]
  }
  base_riskset$eid <- paste0(base_riskset$sender, "_", base_riskset$receiver)
  edgelist$eveid   <- paste0(edgelist$sender, "_", edgelist$receiver)

  if (!is.null(cl)) {
    parallel::clusterExport(
      cl, c("n", "n_nonevents", "edgelist", "riskset_cov2", "base_riskset",
            "covariates", "attach_covariates", ".node_vals", ".find_window"),
      envir = environment()
    )
    pre <- parallel::parLapply(cl, seq(2L, nrow(edgelist)), \(y)
      riskset_cov2(edgelist[seq_len(y - 1L), ], edgelist[y, ],
                   n_nonevents, base_riskset, covariates))
  } else {
    pre <- vector("list", nrow(edgelist))
    for (y in seq(2L, nrow(edgelist)))
      pre[[y]] <- riskset_cov2(edgelist[seq_len(y - 1L), ], edgelist[y, ],
                                n_nonevents, base_riskset, covariates)
  }

  rs <- do.call(rbind.data.frame, pre)
  for (i in which(vapply(rs, is.logical, logical(1L))))
    rs[[i]] <- as.numeric(rs[[i]])

  # NA and non-finite fill applies only to endogenous columns; exogenous
  # covariate columns keep their NAs so the analyst can see missing data.
  exog_cols <- if (!is.null(covariates))
    unlist(lapply(covariates$covariates, function(cv) {
      if (inherits(cv, c("cov_node_static", "cov_node_temporal")) &&
          cv$role == "both")
        c(paste0(cv$name, "_snd"), paste0(cv$name, "_rec"))
      else
        cv$name
    }))
  else
    character(0L)

  endog_cols <- setdiff(names(rs), exog_cols)
  rs[endog_cols] <- lapply(rs[endog_cols], function(x) {
    if (is.numeric(x)) { x[is.na(x)] <- 0; x[!is.finite(x)] <- 0 }
    x
  })
  rs
}
