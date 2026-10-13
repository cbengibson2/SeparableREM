# patch_acc.R
#
# Recompute predictive accuracy for the 205 networks whose acc_sep_init is NA
# in revision_results.Rdata.  The parallel workers silently dropped these via
# the tryCatch in predict_accuracy (model.matrix fails to resolve the formula
# environment in fork-based workers).  Running sequentially here works fine.
#
# Also adds acc_pool_init and acc_pool_cont, which were computed in fit_one_revision
# but never included in the saved res_df columns.

suppressMessages({
  devtools::load_all("llrem", quiet = TRUE)
  library(survival)
})

# ── Formulas ──────────────────────────────────────────────────────────────────

form_init <- Surv(time1, occurred) ~
  RRecInitiation + NOTPSnd + NITPSnd + RecRecencyLog + Imitation + IOviolation +
  strata(time2)

form_cont <- Surv(time1, occurred) ~
  Repetition + NOTPSnd + NITPSnd + RecRecencyLog + Imitation + IOviolation +
  strata(time2)

form_pool <- Surv(time1, occurred) ~
  SndInitiation + RRecInitiation + Repetition + NOTPSnd + NITPSnd +
  RecRecencyLog + Imitation + IOviolation +
  strata(time2)

form_interact <- Surv(time1, occurred) ~
  SndInitiation +
  SndInitiation:RRecInitiation + SndInitiation:Repetition +
  SndInitiation:NOTPSnd + SndInitiation:NITPSnd +
  SndInitiation:RecRecencyLog + SndInitiation:Imitation + SndInitiation:IOviolation +
  RRecInitiation + Repetition + NOTPSnd + NITPSnd +
  RecRecencyLog + Imitation + IOviolation +
  strata(time2)

# ── predict_accuracy ──────────────────────────────────────────────────────────

predict_accuracy <- function(fit, pe_scored, strata_col = "time2", regime = NULL) {
  if (is.null(fit)) return(NA_real_)
  cf <- tryCatch(coef(fit), error = function(e) NULL)
  if (is.null(cf) || all(is.na(cf))) return(NA_real_)
  cf <- cf[!is.na(cf)]
  tt <- tryCatch(
    model.matrix(fit, data = pe_scored),
    error = function(e) NULL
  )
  if (is.null(tt)) return(NA_real_)
  tt <- tt[, names(cf), drop = FALSE]
  lp <- as.numeric(tt %*% cf)
  strata   <- pe_scored[[strata_col]]
  occurred <- pe_scored[["occurred"]]
  snd_init <- if (!is.null(regime)) pe_scored[["SndInitiation"]] else NULL
  unique_strata <- unique(strata)
  correct <- vapply(unique_strata, function(s) {
    idx <- which(strata == s)
    obs <- which(occurred[idx] == 1L)
    if (length(obs) != 1L) return(NA)
    if (!is.null(regime) && snd_init[idx[obs]] != regime) return(NA)
    top <- which.max(lp[idx])
    as.numeric(top == obs)
  }, numeric(1))
  mean(correct, na.rm = TRUE)
}

# ── Load existing results ─────────────────────────────────────────────────────

load("revision_results.Rdata")

needs_patch <- is.na(res_df$acc_sep_init)
paths_to_patch <- rownames(res_df)[needs_patch]
cat(sprintf("Patching %d networks ...\n", sum(needs_patch)))

# Add new columns if absent
if (!"acc_pool_init" %in% names(res_df)) res_df$acc_pool_init <- NA_real_
if (!"acc_pool_cont" %in% names(res_df)) res_df$acc_pool_cont <- NA_real_

# ── Patch loop ────────────────────────────────────────────────────────────────

for (i in seq_along(paths_to_patch)) {
  path <- paths_to_patch[i]
  if (i %% 25 == 0) cat(sprintf("  %d / %d\n", i, length(paths_to_patch)))

  e <- new.env(parent = emptyenv())
  ok <- tryCatch({ load(path, envir = e); TRUE }, error = function(e) FALSE)
  if (!ok) next

  el_raw <- as.data.frame(get("edgelist", envir = e))[, 1:3]
  el     <- tryCatch(prepare_edgelist(el_raw), error = function(e) NULL)
  if (is.null(el) || nrow(el) < 10L) next

  all_nodes <- sort(unique(c(el$sender, el$receiver)))
  base_rs   <- llrem:::get_base_riskset(NULL, all_nodes)
  pe <- tryCatch(
    llrem:::remstats_generate(el, n = length(all_nodes),
                              n_nonevents = 200L, base_riskset = base_rs),
    error = function(e) NULL
  )
  if (is.null(pe)) next
  pe <- llrem:::add_survival_cols(pe)
  if (!"SndInitiation" %in% names(pe)) pe$SndInitiation <- pe$RecInitiation

  pe_init <- pe[pe$SndInitiation == 0, ]
  pe_cont <- pe[pe$SndInitiation == 1, ]

  n_init_ev <- sum(pe_init$occurred)
  n_cont_ev <- sum(pe_cont$occurred)
  if (n_init_ev < 6L || n_cont_ev < 6L) next

  fit_i <- tryCatch(coxph(form_init, data = pe_init), error = function(e) NULL)
  fit_c <- tryCatch(coxph(form_cont, data = pe_cont), error = function(e) NULL)
  fit_p <- tryCatch(coxph(form_pool, data = pe),      error = function(e) NULL)
  fit_x <- tryCatch(coxph(form_interact, data = pe),  error = function(e) NULL)

  acc_i     <- predict_accuracy(fit_i, pe_init)
  acc_c     <- predict_accuracy(fit_c, pe_cont)
  acc_sep   <- if (!is.na(acc_i) && !is.na(acc_c))
    (acc_i * n_init_ev + acc_c * n_cont_ev) / (n_init_ev + n_cont_ev) else NA_real_

  res_df[path, "acc_sep_init"]  <- acc_i
  res_df[path, "acc_sep_cont"]  <- acc_c
  res_df[path, "acc_sep"]       <- acc_sep
  res_df[path, "acc_interact"]  <- predict_accuracy(fit_x, pe)
  res_df[path, "acc_pool_init"] <- predict_accuracy(fit_p, pe, regime = 0L)
  res_df[path, "acc_pool_cont"] <- predict_accuracy(fit_p, pe, regime = 1L)
}

# Also fill acc_pool_init / acc_pool_cont for the 5 networks that already had acc_sep_init
already_ok <- rownames(res_df)[!needs_patch]
cat(sprintf("\nFilling acc_pool_init/cont for %d already-ok networks ...\n", length(already_ok)))

for (path in already_ok) {
  e <- new.env(parent = emptyenv())
  ok <- tryCatch({ load(path, envir = e); TRUE }, error = function(e) FALSE)
  if (!ok) next
  el_raw <- as.data.frame(get("edgelist", envir = e))[, 1:3]
  el     <- tryCatch(prepare_edgelist(el_raw), error = function(e) NULL)
  if (is.null(el) || nrow(el) < 10L) next
  all_nodes <- sort(unique(c(el$sender, el$receiver)))
  base_rs   <- llrem:::get_base_riskset(NULL, all_nodes)
  pe <- tryCatch(
    llrem:::remstats_generate(el, n = length(all_nodes),
                              n_nonevents = 200L, base_riskset = base_rs),
    error = function(e) NULL
  )
  if (is.null(pe)) next
  pe <- llrem:::add_survival_cols(pe)
  if (!"SndInitiation" %in% names(pe)) pe$SndInitiation <- pe$RecInitiation
  fit_p <- tryCatch(coxph(form_pool, data = pe), error = function(e) NULL)
  res_df[path, "acc_pool_init"] <- predict_accuracy(fit_p, pe, regime = 0L)
  res_df[path, "acc_pool_cont"] <- predict_accuracy(fit_p, pe, regime = 1L)
}

# ── Summary ───────────────────────────────────────────────────────────────────

cat(sprintf("\n=== After patch ===\n"))
cat(sprintf("acc_sep_init  non-NA: %d / %d\n", sum(!is.na(res_df$acc_sep_init)),  nrow(res_df)))
cat(sprintf("acc_sep_cont  non-NA: %d / %d\n", sum(!is.na(res_df$acc_sep_cont)),  nrow(res_df)))
cat(sprintf("acc_sep       non-NA: %d / %d\n", sum(!is.na(res_df$acc_sep)),       nrow(res_df)))
cat(sprintf("acc_interact  non-NA: %d / %d\n", sum(!is.na(res_df$acc_interact)),  nrow(res_df)))
cat(sprintf("acc_pool_init non-NA: %d / %d\n", sum(!is.na(res_df$acc_pool_init)), nrow(res_df)))
cat(sprintf("acc_pool_cont non-NA: %d / %d\n", sum(!is.na(res_df$acc_pool_cont)), nrow(res_df)))

cat(sprintf("\nMean acc_sep_init:  %.4f\n", mean(res_df$acc_sep_init,  na.rm = TRUE)))
cat(sprintf("Mean acc_sep_cont:  %.4f\n", mean(res_df$acc_sep_cont,  na.rm = TRUE)))
cat(sprintf("Mean acc_pool:      %.4f\n", mean(res_df$acc_pool,      na.rm = TRUE)))
cat(sprintf("Mean acc_pool_init: %.4f\n", mean(res_df$acc_pool_init, na.rm = TRUE)))
cat(sprintf("Mean acc_pool_cont: %.4f\n", mean(res_df$acc_pool_cont, na.rm = TRUE)))

# ── Save ──────────────────────────────────────────────────────────────────────

save(res_df, file = "revision_results.Rdata")
cat("\nSaved to revision_results.Rdata\n")
