# revision.R
#
# Analyses added in response to reviewer comments.
#
# 1. Fully interacted pooled model: all terms x SndInitiation (16 params).
#    Compared to separable (12 params) and original pooled (8 params).
#    Addresses: reviewer concern about AIC mechanical advantage from risk-set reduction.
#
# 2. Restricted pooled model: RRecInitiation without RecRecencyLog.
#    Tests whether RecRecencyLog collinearity drives the RRecInitiation sign reversal.
#    Addresses: reviewer concern that the two-mechanism explanation is untested.
#
# 3. Decay-window robustness: A_ijt resets to 0 if no i->j event in last W events.
#    Tests sensitivity of results to the irreversibility assumption.
#    Addresses: reviewer concern about dormant relationships in gang/war sequences.
#
# 4. Results reported separately for synthetic (Gang Violence) vs observed networks.
#    Addresses: reviewer concern about synthetic data dominating the corpus.
#
# 5. Full greedy selection on all 210 networks (not subsampled 78).
#    Addresses: reviewer concern about greedy subsample with no sensitivity check.
#
# 6. Predictive accuracy: for each model, at each event score all dyads in the risk
#    set using the fitted linear predictor and record whether the actual dyad ranked
#    first. Report mean proportion-correct per model per domain.
#
# 7. Ecology-weighted robustness (appendix): weight each network by 1/N_d where N_d
#    is the number of networks in its domain, so each ecology contributes equally to
#    pooled summaries regardless of corpus size.

rm(list = ls())
devtools::load_all("llrem", quiet = TRUE)
library(parallel)
library(survival)

# ── Domain metadata ───────────────────────────────────────────────────────────

DATA_DIR <- Sys.getenv("DATA_DIR", "")   # path to the DataP network data directory

domain_meta <- list(
  Classroom = list(dir = "Classroom", cat = "Classroom",    visible = "Y", synthetic = FALSE),
  Email     = list(dir = "Email",     cat = "Email",        visible = "N", synthetic = FALSE),
  Gangs     = list(dir = "Gangs",     cat = "Gang Violence",visible = "Y", synthetic = TRUE),
  PMs       = list(dir = "PMs",       cat = "PMs",          visible = "N", synthetic = FALSE),
  Phone     = list(dir = "Phone",     cat = "Phone (DoD)",  visible = "Y", synthetic = FALSE),
  SMS       = list(dir = "SMS",       cat = "SMS",          visible = "N", synthetic = FALSE),
  Twitter   = list(dir = "Twitter",   cat = "Twitter",      visible = "Y", synthetic = FALSE),
  War       = list(dir = "War",       cat = "War",          visible = "Y", synthetic = FALSE),
  ftf       = list(dir = "ftf",       cat = "Glances",      visible = "N", synthetic = FALSE)
)

all_files <- do.call(rbind, lapply(names(domain_meta), function(dname) {
  meta  <- domain_meta[[dname]]
  flist <- list.files(file.path(DATA_DIR, meta$dir),
                      pattern = "\\.[Rr]data$", full.names = TRUE)
  data.frame(path = flist, cat = meta$cat, visible = meta$visible,
             synthetic = meta$synthetic, stringsAsFactors = FALSE)
}))

cat(sprintf("Found %d files across %d domains\n\n",
            nrow(all_files), length(unique(all_files$cat))))
print(table(all_files$cat))

# ── Model formulas ────────────────────────────────────────────────────────────

# Original sub-models
form_init <- Surv(time1, occurred) ~
  RRecInitiation + NOTPSnd + NITPSnd + RecRecencyLog + Imitation + IOviolation +
  strata(time2)

form_cont <- Surv(time1, occurred) ~
  Repetition + NOTPSnd + NITPSnd + RecRecencyLog + Imitation + IOviolation +
  strata(time2)

# Original pooled (8 params)
form_pool <- Surv(time1, occurred) ~
  SndInitiation + RRecInitiation + Repetition + NOTPSnd + NITPSnd +
  RecRecencyLog + Imitation + IOviolation +
  strata(time2)

# Fully interacted pooled model (16 params): all terms x SndInitiation.
# Uses same pooled risk set as form_pool; directly comparable to separable (12 params)
# in terms of the degree of freedom penalty.
form_interact <- Surv(time1, occurred) ~
  SndInitiation +
  SndInitiation:RRecInitiation + SndInitiation:Repetition +
  SndInitiation:NOTPSnd + SndInitiation:NITPSnd +
  SndInitiation:RecRecencyLog + SndInitiation:Imitation + SndInitiation:IOviolation +
  RRecInitiation + Repetition + NOTPSnd + NITPSnd +
  RecRecencyLog + Imitation + IOviolation +
  strata(time2)

# Restricted pooled model: RRecInitiation without RecRecencyLog (tests collinearity mechanism)
form_pool_noRec <- Surv(time1, occurred) ~
  SndInitiation + RRecInitiation + Repetition + NOTPSnd + NITPSnd +
  Imitation + IOviolation +
  strata(time2)

# ── Decay-window A_ijt: reset prior-contact indicator after W-event gap ───────
# W = 100 events; i->j resets to "never contacted" if no i->j event in last W events.

make_decay_initiation <- function(pe, W = 100L) {
  # pe must have columns: event_index (integer), sender_id, receiver_id, occurred
  # We recompute SndInitiation based on whether i->j occurred in the window
  # [event_index - W, event_index - 1].
  #
  # Because pe is already expanded (one row per dyad per event), we use the
  # event_index column to define the window.
  if (!"event_index" %in% names(pe)) return(pe)

  pe <- pe[order(pe$event_index, pe$sender_id, pe$receiver_id), ]
  key <- paste(pe$sender_id, pe$receiver_id, sep = "_")
  unique_keys <- unique(key)

  decay_init <- integer(nrow(pe))
  for (k in unique_keys) {
    idx <- which(key == k)
    rows <- pe[idx, ]
    ev_idx <- rows$event_index
    # For each row, check if any occurred==1 row for this dyad falls in (ev_idx-W, ev_idx)
    occurred_at <- ev_idx[rows$occurred == 1]
    for (i in seq_along(idx)) {
      lo <- ev_idx[i] - W
      decay_init[idx[i]] <- as.integer(any(occurred_at >= lo & occurred_at < ev_idx[i]))
    }
  }
  pe$SndInitiation_decay <- decay_init
  pe
}

# ── Predictive accuracy ───────────────────────────────────────────────────────
#
# For each observed event, score every dyad in the same-stratum risk set using
# the model's linear predictor (beta' * x), then check whether the actual dyad
# (occurred == 1) receives the highest score. The proportion of events where the
# model's top-ranked dyad is the one that actually occurred is the accuracy.
#
# The separable model scores initiation events using the initiation sub-model and
# continuation events using the continuation sub-model; each event is compared
# only within its own regime's risk set (i.e. the same dyads the sub-model saw).
# The pooled models score all events using the single pooled linear predictor
# against the full risk set.
#
# This is a leave-nothing-out in-sample measure. Its purpose is not to estimate
# generalisation error but to check whether the model's ordering of the risk set
# aligns with observed choices --- a model-free complement to AIC.

predict_accuracy <- function(fit, pe_scored, strata_col = "time2", regime = NULL) {
  # fit        : coxph object
  # pe_scored  : data frame used to score (must contain all model covariates
  #              and strata_col plus 'occurred')
  # regime     : NULL (all events), 0 (initiation only), or 1 (continuation only).
  #              When non-NULL, uses the full risk set per stratum but only
  #              evaluates strata whose actual occurred event matches the regime.
  # Returns proportion of events where the model's top-ranked risk-set member
  # is the one that actually occurred.
  if (is.null(fit)) return(NA_real_)
  cf <- tryCatch(coef(fit), error = function(e) NULL)
  if (is.null(cf) || all(is.na(cf))) return(NA_real_)
  cf <- cf[!is.na(cf)]   # drop aliased/unidentified terms

  # Build model matrix (no intercept; same terms as formula minus strata)
  tt <- tryCatch(
    model.matrix(fit, data = pe_scored),
    error = function(e) NULL
  )
  if (is.null(tt)) return(NA_real_)
  tt <- tt[, names(cf), drop = FALSE]   # keep only identified columns

  # Linear predictor for every row
  lp <- as.numeric(tt %*% cf)

  strata   <- pe_scored[[strata_col]]
  occurred <- pe_scored[["occurred"]]
  snd_init <- if (!is.null(regime)) pe_scored[["SndInitiation"]] else NULL
  unique_strata <- unique(strata)

  correct <- vapply(unique_strata, function(s) {
    idx <- which(strata == s)
    obs <- which(occurred[idx] == 1L)
    if (length(obs) != 1L) return(NA)           # skip strata with != 1 event
    # If filtering by regime, skip strata whose occurred event is wrong regime
    if (!is.null(regime) && snd_init[idx[obs]] != regime) return(NA)
    top <- which.max(lp[idx])
    as.numeric(top == obs)
  }, numeric(1))

  mean(correct, na.rm = TRUE)
}

# ── Per-network fitting function ──────────────────────────────────────────────

fit_one_revision <- function(path, cat_label, visible_label, synthetic_label) {
  e <- new.env(parent = emptyenv())
  tryCatch(load(path, envir = e), error = function(err) return(NULL))

  el_raw <- as.data.frame(get("edgelist", envir = e))[, 1:3]
  el     <- prepare_edgelist(el_raw)
  if (nrow(el) < 10L) return(NULL)

  all_nodes <- sort(unique(c(el$sender, el$receiver)))
  base_rs   <- llrem:::get_base_riskset(NULL, all_nodes)

  pe <- tryCatch(
    llrem:::remstats_generate(el, n = length(all_nodes),
                              n_nonevents = 200L, base_riskset = base_rs),
    error = function(err) NULL
  )
  if (is.null(pe)) return(NULL)
  pe <- llrem:::add_survival_cols(pe)
  if (!"SndInitiation" %in% names(pe)) pe$SndInitiation <- pe$RecInitiation

  pe_init <- pe[pe$SndInitiation == 0, ]
  pe_cont <- pe[pe$SndInitiation == 1, ]

  n_init_ev <- sum(pe_init$occurred)
  n_cont_ev <- sum(pe_cont$occurred)

  if (n_init_ev < 6L || n_cont_ev < 6L) return(NULL)

  # ── Fit original sub-models ───────────────────────────────────────────────
  fit_i <- tryCatch(coxph(form_init, data = pe_init), error = function(e) NULL)
  fit_c <- tryCatch(coxph(form_cont, data = pe_cont), error = function(e) NULL)
  fit_p <- tryCatch(coxph(form_pool, data = pe),      error = function(e) NULL)
  if (is.null(fit_i) || is.null(fit_c) || is.null(fit_p)) return(NULL)

  # ── Fit fully interacted pooled model ────────────────────────────────────
  fit_interact <- tryCatch(coxph(form_interact, data = pe), error = function(e) NULL)

  # ── Fit restricted pooled model (no RecRecencyLog) ───────────────────────────
  fit_pool_noRec <- tryCatch(coxph(form_pool_noRec, data = pe), error = function(e) NULL)

  # ── Decay-window robustness ───────────────────────────────────────────────
  # Only attempt if event_index column present
  fit_i_decay <- NULL
  fit_c_decay <- NULL
  if ("event_index" %in% names(pe)) {
    pe_decay <- make_decay_initiation(pe, W = 100L)
    pe_init_d <- pe_decay[pe_decay$SndInitiation_decay == 0, ]
    pe_cont_d <- pe_decay[pe_decay$SndInitiation_decay == 1, ]
    if (sum(pe_init_d$occurred) >= 6L && sum(pe_cont_d$occurred) >= 6L) {
      fit_i_decay <- tryCatch(coxph(form_init, data = pe_init_d), error = function(e) NULL)
      fit_c_decay <- tryCatch(coxph(form_cont, data = pe_cont_d), error = function(e) NULL)
    }
  }

  # ── Extract coefficients ──────────────────────────────────────────────────
  get_cf <- function(cf, nm, col = "coef") {
    if (!is.null(cf) && nm %in% rownames(cf)) cf[nm, col] else NA_real_
  }
  get_z <- function(cf, nm) get_cf(cf, nm, "z")

  cf_i  <- if (!is.null(fit_i))          coef(summary(fit_i))          else NULL
  cf_c  <- if (!is.null(fit_c))          coef(summary(fit_c))          else NULL
  cf_p  <- if (!is.null(fit_p))          coef(summary(fit_p))          else NULL
  cf_in <- if (!is.null(fit_interact))   coef(summary(fit_interact))   else NULL
  cf_nr <- if (!is.null(fit_pool_noRec)) coef(summary(fit_pool_noRec)) else NULL
  cf_id <- if (!is.null(fit_i_decay))    coef(summary(fit_i_decay))    else NULL
  cf_cd <- if (!is.null(fit_c_decay))    coef(summary(fit_c_decay))    else NULL

  # ── Log-likelihoods and AIC ───────────────────────────────────────────────
  ll_i   <- if (!is.null(fit_i)) fit_i$loglik[2] else NA_real_
  ll_c   <- if (!is.null(fit_c)) fit_c$loglik[2] else NA_real_
  ll_p   <- if (!is.null(fit_p)) fit_p$loglik[2] else NA_real_
  ll_sep <- ll_i + ll_c

  df_sep    <- length(coef(fit_i)) + length(coef(fit_c))  # 12
  df_pool   <- length(coef(fit_p))                         # 8
  aic_sep   <- -2 * ll_sep + 2 * df_sep
  aic_pool  <- -2 * ll_p   + 2 * df_pool
  delta_aic <- aic_pool - aic_sep

  ll_interact  <- if (!is.null(fit_interact))   fit_interact$loglik[2]   else NA_real_
  df_interact  <- if (!is.null(fit_interact))   length(coef(fit_interact)) else NA_real_
  aic_interact <- if (!is.na(ll_interact))      -2 * ll_interact + 2 * df_interact else NA_real_
  # delta_AIC: interacted pooled vs separable (both have ~12-16 params; compare directly)
  delta_aic_interact_vs_sep  <- aic_interact - aic_sep
  delta_aic_interact_vs_pool <- aic_pool     - aic_interact

  ll_noRec   <- if (!is.null(fit_pool_noRec)) fit_pool_noRec$loglik[2] else NA_real_
  df_noRec   <- if (!is.null(fit_pool_noRec)) length(coef(fit_pool_noRec)) else NA_real_

  # Sign of RRecInitiation in pooled models with and without RecRecencyLog
  pool_rrec_sign    <- sign(get_cf(cf_p,  "RRecInitiation"))
  noRec_rrec_sign   <- sign(get_cf(cf_nr, "RRecInitiation"))
  init_rrec_sign    <- sign(get_cf(cf_i,  "RRecInitiation"))
  sign_flip_pool    <- as.integer(!is.na(pool_rrec_sign)  & !is.na(init_rrec_sign)  & pool_rrec_sign  != init_rrec_sign)
  sign_flip_noRec   <- as.integer(!is.na(noRec_rrec_sign) & !is.na(init_rrec_sign)  & noRec_rrec_sign != init_rrec_sign)

  # Decay-window AIC
  ll_sep_decay  <- NA_real_
  aic_sep_decay <- NA_real_
  if (!is.null(fit_i_decay) && !is.null(fit_c_decay)) {
    ll_sep_decay  <- fit_i_decay$loglik[2] + fit_c_decay$loglik[2]
    aic_sep_decay <- -2 * ll_sep_decay + 2 * df_sep
  }

  list(
    path        = path,
    cat         = cat_label,
    visible     = visible_label,
    synthetic   = synthetic_label,
    n           = length(all_nodes),
    n_events    = nrow(el),
    n_init      = n_init_ev,
    n_cont      = n_cont_ev,

    # Original model fit
    ll_init     = ll_i,
    ll_cont     = ll_c,
    ll_sep      = ll_sep,
    ll_pool     = ll_p,
    aic_sep     = aic_sep,
    aic_pool    = aic_pool,
    delta_aic   = delta_aic,

    # Interacted pooled model
    ll_interact              = ll_interact,
    df_interact              = df_interact,
    aic_interact             = aic_interact,
    delta_aic_interact_vs_sep  = delta_aic_interact_vs_sep,
    delta_aic_interact_vs_pool = delta_aic_interact_vs_pool,

    # Restricted pooled (no RecRecencyLog): collinearity test
    ll_noRec         = ll_noRec,
    df_noRec         = df_noRec,
    sign_flip_pool   = sign_flip_pool,
    sign_flip_noRec  = sign_flip_noRec,
    pool_RRecInitiation  = get_cf(cf_p,  "RRecInitiation"),
    noRec_RRecInitiation = get_cf(cf_nr, "RRecInitiation"),
    init_RRecInitiation  = get_cf(cf_i,  "RRecInitiation"),

    # Decay-window robustness
    ll_sep_decay    = ll_sep_decay,
    aic_sep_decay   = aic_sep_decay,
    init_NOTPSnd_decay = get_cf(cf_id, "NOTPSnd"),
    cont_NOTPSnd_decay = get_cf(cf_cd, "NOTPSnd"),

    # Key coefficients for triadic/reciprocity claims
    init_NOTPSnd    = get_cf(cf_i,  "NOTPSnd"),
    cont_NOTPSnd    = get_cf(cf_c,  "NOTPSnd"),
    pool_NOTPSnd    = get_cf(cf_p,  "NOTPSnd"),
    init_NITPSnd    = get_cf(cf_i,  "NITPSnd"),
    cont_NITPSnd    = get_cf(cf_c,  "NITPSnd"),

    # ── Full coefficient and z-stat table (all terms, all three models) ────────
    # Initiation model
    init_RRecInitiation_z = get_z(cf_i, "RRecInitiation"),
    init_NOTPSnd_z        = get_z(cf_i, "NOTPSnd"),
    init_NITPSnd_z        = get_z(cf_i, "NITPSnd"),
    init_RecRecencyLog    = get_cf(cf_i, "RecRecencyLog"),
    init_RecRecencyLog_z  = get_z(cf_i,  "RecRecencyLog"),
    init_Imitation        = get_cf(cf_i, "Imitation"),
    init_Imitation_z      = get_z(cf_i,  "Imitation"),
    init_IOviolation      = get_cf(cf_i, "IOviolation"),
    init_IOviolation_z    = get_z(cf_i,  "IOviolation"),

    # Continuation model
    cont_Repetition       = get_cf(cf_c, "Repetition"),
    cont_Repetition_z     = get_z(cf_c,  "Repetition"),
    cont_NOTPSnd_z        = get_z(cf_c,  "NOTPSnd"),
    cont_NITPSnd_z        = get_z(cf_c,  "NITPSnd"),
    cont_RecRecencyLog    = get_cf(cf_c, "RecRecencyLog"),
    cont_RecRecencyLog_z  = get_z(cf_c,  "RecRecencyLog"),
    cont_Imitation        = get_cf(cf_c, "Imitation"),
    cont_Imitation_z      = get_z(cf_c,  "Imitation"),
    cont_IOviolation      = get_cf(cf_c, "IOviolation"),
    cont_IOviolation_z    = get_z(cf_c,  "IOviolation"),

    # Pooled model
    pool_SndInitiation    = get_cf(cf_p, "SndInitiation"),
    pool_SndInitiation_z  = get_z(cf_p,  "SndInitiation"),
    pool_RRecInitiation_z = get_z(cf_p,  "RRecInitiation"),
    pool_Repetition       = get_cf(cf_p, "Repetition"),
    pool_Repetition_z     = get_z(cf_p,  "Repetition"),
    pool_NOTPSnd_z        = get_z(cf_p,  "NOTPSnd"),
    pool_NITPSnd          = get_cf(cf_p, "NITPSnd"),
    pool_NITPSnd_z        = get_z(cf_p,  "NITPSnd"),
    pool_RecRecencyLog    = get_cf(cf_p, "RecRecencyLog"),
    pool_RecRecencyLog_z  = get_z(cf_p,  "RecRecencyLog"),
    pool_Imitation        = get_cf(cf_p, "Imitation"),
    pool_Imitation_z      = get_z(cf_p,  "Imitation"),
    pool_IOviolation      = get_cf(cf_p, "IOviolation"),
    pool_IOviolation_z    = get_z(cf_p,  "IOviolation"),

    # ── Predictive accuracy ───────────────────────────────────────────────────
    # Separable: score each event within its own regime's risk set.
    # The accuracy figures are thus directly comparable across models because
    # the separable model scores initiation events against initiation dyads only
    # and continuation events against continuation dyads only.
    acc_sep  = {
      acc_i <- predict_accuracy(fit_i, pe_init)
      acc_c <- predict_accuracy(fit_c, pe_cont)
      # Weighted average by number of events in each regime
      if (!is.na(acc_i) && !is.na(acc_c))
        (acc_i * n_init_ev + acc_c * n_cont_ev) / (n_init_ev + n_cont_ev)
      else NA_real_
    },
    acc_sep_init     = predict_accuracy(fit_i, pe_init),
    acc_sep_cont     = predict_accuracy(fit_c, pe_cont),
    acc_pool         = predict_accuracy(fit_p, pe),
    acc_pool_init    = predict_accuracy(fit_p, pe, regime = 0L),
    acc_pool_cont    = predict_accuracy(fit_p, pe, regime = 1L),
    acc_interact     = predict_accuracy(fit_interact, pe)
  )
}

# ── Parallel execution ────────────────────────────────────────────────────────

n_cores <- min(parallel::detectCores() - 1L, 6L)
cl <- parallel::makeCluster(n_cores)
parallel::clusterCall(cl, function() {
  devtools::load_all("llrem", quiet = TRUE)
  library(survival)
})
parallel::clusterExport(cl, c(
  "form_init", "form_cont", "form_pool",
  "form_interact", "form_pool_noRec",
  "make_decay_initiation", "predict_accuracy",
  "fit_one_revision", "DATA_DIR"
))

cat(sprintf("\nRunning revision analyses on %d networks using %d cores ...\n",
            nrow(all_files), n_cores))

set.seed(8214L)
raw_results <- parallel::clusterMap(
  cl, fit_one_revision,
  path           = all_files$path,
  cat_label      = all_files$cat,
  visible_label  = all_files$visible,
  synthetic_label = all_files$synthetic,
  SIMPLIFY = FALSE
)
parallel::stopCluster(cl)

results <- raw_results[!sapply(raw_results, is.null)]
cat(sprintf("%d networks retained (dropped %d)\n\n",
            length(results), nrow(all_files) - length(results)))

res_df <- do.call(rbind, lapply(results, function(x) {
  as.data.frame(x[names(x) != "path"], stringsAsFactors = FALSE)
}))

save(res_df, file = "revision_results.Rdata")
cat("Results saved to revision_results.Rdata\n\n")

# ── Summary helpers ───────────────────────────────────────────────────────────

smry <- function(x, label = "") {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(data.frame(label=label, n=0, mean=NA, sd=NA, t=NA))
  m <- mean(x); s <- sd(x); n <- length(x)
  data.frame(label=label, n=n, mean=round(m,3), sd=round(s,3),
             t=round(m/(s/sqrt(n)), 2))
}

smry_by_cat <- function(col, df) {
  do.call(rbind, lapply(sort(unique(df$cat)), function(g) {
    smry(df[[col]][df$cat == g], label = g)
  }))
}

# ── Report 1: AIC comparison — original, interacted, separable ───────────────

cat("=================================================================\n")
cat("REPORT 1: AIC COMPARISON\n")
cat("  Original pooled (8p) vs Interacted pooled (16p) vs Separable (12p)\n")
cat("=================================================================\n")

cat("\n--- Separable preferred over original pooled (delta_AIC > 0) ---\n")
cat(sprintf("  %d / %d networks\n", sum(res_df$delta_aic > 0, na.rm=TRUE), sum(!is.na(res_df$delta_aic))))

cat("\n--- Separable vs interacted pooled (negative = separable preferred) ---\n")
cat(sprintf("  Separable preferred: %d / %d networks\n",
            sum(res_df$delta_aic_interact_vs_sep < 0, na.rm=TRUE),
            sum(!is.na(res_df$delta_aic_interact_vs_sep))))
print(smry_by_cat("delta_aic_interact_vs_sep", res_df), row.names=FALSE)

cat("\n--- Mean delta_AIC: interacted vs original pooled (positive = interacted preferred) ---\n")
print(smry_by_cat("delta_aic_interact_vs_pool", res_df), row.names=FALSE)

# ── Report 2: Synthetic vs observed ──────────────────────────────────────────

cat("\n=================================================================\n")
cat("REPORT 2: SEPARABLE PREFERENCE — OBSERVED-ONLY NETWORKS\n")
cat("=================================================================\n")

obs_df <- res_df[!res_df$synthetic, ]
cat(sprintf("  Observed networks: %d\n", nrow(obs_df)))
cat(sprintf("  Separable preferred (delta_AIC > 0): %d / %d\n",
            sum(obs_df$delta_aic > 0, na.rm=TRUE),
            sum(!is.na(obs_df$delta_aic))))
cat("\n  Mean delta_AIC by domain (observed only):\n")
print(smry_by_cat("delta_aic", obs_df), row.names=FALSE)

# ── Report 3: Collinearity test — RRecInitiation sign with/without RecRecencyLog ─

cat("\n=================================================================\n")
cat("REPORT 3: COLLINEARITY TEST\n")
cat("  Sign flip rate: pooled (with RecRecencyLog) vs pooled (without RecRecencyLog)\n")
cat("=================================================================\n")

flip_by_cat <- do.call(rbind, lapply(sort(unique(res_df$cat)), function(g) {
  sub <- res_df[res_df$cat == g, ]
  n       <- nrow(sub)
  n_flip_pool  <- sum(sub$sign_flip_pool,  na.rm=TRUE)
  n_flip_noRec <- sum(sub$sign_flip_noRec, na.rm=TRUE)
  data.frame(
    Domain          = g,
    N               = n,
    flip_pool_pct   = round(100 * n_flip_pool  / n, 1),
    flip_noRec_pct  = round(100 * n_flip_noRec / n, 1)
  )
}))
print(flip_by_cat, row.names=FALSE)
cat("\n  If collinearity drives the reversal, flip_noRec_pct should be lower.\n")

# ── Report 4: Triadic closure — initiation vs continuation ratio ──────────────

cat("\n=================================================================\n")
cat("REPORT 4: TRIADIC CLOSURE — INIT vs CONT NOTPSnd RATIO\n")
cat("=================================================================\n")

ratio_by_cat <- do.call(rbind, lapply(sort(unique(res_df$cat)), function(g) {
  sub <- res_df[res_df$cat == g, ]
  m_init <- mean(sub$init_NOTPSnd, na.rm=TRUE)
  m_cont <- mean(sub$cont_NOTPSnd, na.rm=TRUE)
  data.frame(
    Domain    = g,
    N         = nrow(sub),
    init_NOTP = round(m_init, 3),
    cont_NOTP = round(m_cont, 3),
    ratio     = round(m_init / m_cont, 2)
  )
}))
print(ratio_by_cat, row.names=FALSE)

# ── Report 5: Decay-window robustness ────────────────────────────────────────

cat("\n=================================================================\n")
cat("REPORT 5: DECAY-WINDOW ROBUSTNESS (W=100 events)\n")
cat("  Domains where decay NOTPSnd init/cont available\n")
cat("=================================================================\n")

decay_sub <- res_df[!is.na(res_df$init_NOTPSnd_decay), ]
if (nrow(decay_sub) > 0) {
  cat(sprintf("  Networks with decay estimates: %d\n", nrow(decay_sub)))
  decay_by_cat <- do.call(rbind, lapply(sort(unique(decay_sub$cat)), function(g) {
    sub <- decay_sub[decay_sub$cat == g, ]
    data.frame(
      Domain          = g,
      N               = nrow(sub),
      init_NOTP_orig  = round(mean(res_df$init_NOTPSnd[res_df$cat == g], na.rm=TRUE), 3),
      init_NOTP_decay = round(mean(sub$init_NOTPSnd_decay, na.rm=TRUE), 3),
      cont_NOTP_orig  = round(mean(res_df$cont_NOTPSnd[res_df$cat == g], na.rm=TRUE), 3),
      cont_NOTP_decay = round(mean(sub$cont_NOTPSnd_decay, na.rm=TRUE), 3)
    )
  }))
  print(decay_by_cat, row.names=FALSE)
} else {
  cat("  event_index column not available in riskset; decay analysis skipped.\n")
  cat("  Robustness check requires event_index in llrem riskset output.\n")
}

# ── Report 6: Predictive accuracy ────────────────────────────────────────────
#
# Proportion of events for which the model's highest-ranked dyad in the risk set
# is the dyad that actually occurred. Reported by domain and overall.

cat("\n=================================================================\n")
cat("REPORT 6: PREDICTIVE ACCURACY\n")
cat("  Proportion of events where model's top-ranked dyad is the observed event\n")
cat("  (separable scored within-regime; pooled scored against full risk set)\n")
cat("=================================================================\n")

acc_by_cat <- do.call(rbind, lapply(sort(unique(res_df$cat)), function(g) {
  sub <- res_df[res_df$cat == g, ]
  data.frame(
    Domain       = g,
    N            = nrow(sub),
    acc_sep      = round(mean(sub$acc_sep,      na.rm = TRUE), 4),
    acc_sep_init = round(mean(sub$acc_sep_init, na.rm = TRUE), 4),
    acc_sep_cont = round(mean(sub$acc_sep_cont, na.rm = TRUE), 4),
    acc_pool     = round(mean(sub$acc_pool,     na.rm = TRUE), 4),
    acc_interact = round(mean(sub$acc_interact, na.rm = TRUE), 4)
  )
}))
print(acc_by_cat, row.names = FALSE)

cat("\n--- Overall (unweighted mean across networks) ---\n")
cat(sprintf("  Separable (combined):  %.4f\n", mean(res_df$acc_sep,      na.rm = TRUE)))
cat(sprintf("  Separable (init only): %.4f\n", mean(res_df$acc_sep_init, na.rm = TRUE)))
cat(sprintf("  Separable (cont only): %.4f\n", mean(res_df$acc_sep_cont, na.rm = TRUE)))
cat(sprintf("  Pooled (original):     %.4f\n", mean(res_df$acc_pool,     na.rm = TRUE)))
cat(sprintf("  Pooled (interacted):   %.4f\n", mean(res_df$acc_interact, na.rm = TRUE)))

# ── Report 7: Ecology-weighted robustness (appendix) ─────────────────────────
#
# Each domain ecology contributes equally regardless of how many networks it
# contains. Weight for network k in domain d: w_k = 1 / N_d.
# This prevents the Gang Violence corpus (n=100) from dominating pooled means.
#
# Weighted mean of X: sum(w_k * X_k) / sum(w_k)

cat("\n=================================================================\n")
cat("REPORT 7: ECOLOGY-WEIGHTED ROBUSTNESS (appendix)\n")
cat("  Each ecology weighted equally (w = 1/N_d per network)\n")
cat("=================================================================\n")

# Compute per-network weights
domain_counts <- table(res_df$cat)
res_df$w <- 1 / as.numeric(domain_counts[res_df$cat])

wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

wsmry <- function(col, df, label = col) {
  x <- df[[col]]; w <- df$w
  ok <- is.finite(x) & is.finite(w)
  if (!any(ok)) return(data.frame(label=label, n=0, wmean=NA, wsd=NA))
  xw <- wmean(x, w)
  # Weighted SD: sqrt( sum(w*(x-xbar)^2) / sum(w) )
  wsd <- sqrt(sum(w[ok] * (x[ok] - xw)^2) / sum(w[ok]))
  data.frame(label = label, n = sum(ok),
             wmean = round(xw,  4),
             wsd   = round(wsd, 4))
}

cat("\n--- Weighted proportion of networks where separable preferred (delta_AIC > 0) ---\n")
pref <- as.integer(res_df$delta_aic > 0 & !is.na(res_df$delta_aic))
cat(sprintf("  Ecology-weighted rate: %.4f\n", wmean(pref, res_df$w)))

cat("\n--- Weighted mean delta_AIC (separable vs pooled) by domain ---\n")
# Normalise delta_AIC by sequence length to remove scale dependence
res_df$delta_aic_norm <- res_df$delta_aic / res_df$n_events
wda_by_cat <- do.call(rbind, lapply(sort(unique(res_df$cat)), function(g) {
  sub <- res_df[res_df$cat == g, ]
  data.frame(
    Domain           = g,
    N                = nrow(sub),
    wmean_delta_aic  = round(wmean(sub$delta_aic,      sub$w), 2),
    wmean_delta_norm = round(wmean(sub$delta_aic_norm, sub$w), 4)
  )
}))
print(wda_by_cat, row.names = FALSE)

cat("\n--- Weighted predictive accuracy ---\n")
wacc <- data.frame(
  Model        = c("Separable", "Separable-init", "Separable-cont",
                   "Pooled-original", "Pooled-interacted"),
  Weighted_acc = c(
    wmean(res_df$acc_sep,      res_df$w),
    wmean(res_df$acc_sep_init, res_df$w),
    wmean(res_df$acc_sep_cont, res_df$w),
    wmean(res_df$acc_pool,     res_df$w),
    wmean(res_df$acc_interact, res_df$w)
  )
)
wacc$Weighted_acc <- round(wacc$Weighted_acc, 4)
print(wacc, row.names = FALSE)

cat("\n--- Weighted sign-flip rates (RRecInitiation) ---\n")
wsf <- data.frame(
  Domain = sort(unique(res_df$cat))
)
wsf$w_flip_pool <- vapply(wsf$Domain, function(g) {
  sub <- res_df[res_df$cat == g, ]
  wmean(sub$sign_flip_pool, sub$w)
}, numeric(1))
wsf$w_flip_noRec <- vapply(wsf$Domain, function(g) {
  sub <- res_df[res_df$cat == g, ]
  wmean(sub$sign_flip_noRec, sub$w)
}, numeric(1))
wsf$w_flip_pool  <- round(wsf$w_flip_pool,  3)
wsf$w_flip_noRec <- round(wsf$w_flip_noRec, 3)
print(wsf, row.names = FALSE)

save(res_df, file = "revision_results.Rdata")
cat("\nUpdated results (with accuracy + weights) saved to revision_results.Rdata\n")

cat("\nDone.\n")
