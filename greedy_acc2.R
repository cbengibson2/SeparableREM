# greedy_acc2.R
#
# Faster companion to greedy_acc.R.
#
# For the SEPARABLE model: reads the previously stored greedy_sep_results.Rdata
# and reconstructs each network's final selected model by applying the stopping
# rule (delta_aic < 0) to the recorded path. Only that final model is refit;
# no forward search is repeated.
#
# For the POOLED model: no prior greedy results exist, so greedy forward
# selection is run fresh (same algorithm and stopping rule as greedy_acc.R).
#
# Saves coefficients, SEs, z-statistics, and predictive accuracy to
# greedy_acc2_results.Rdata.

rm(list = ls())
devtools::load_all("llrem", quiet = TRUE)
library(parallel)
library(survival)

# ── Paths ─────────────────────────────────────────────────────────────────────

DATA_DIR    <- Sys.getenv("DATA_DIR",    "")   # path to the DataP network data directory
OUT_DIR     <- Sys.getenv("OUT_DIR",     ".")  # output directory (default: project root)
GREEDY_FILE <- Sys.getenv("GREEDY_FILE",
                file.path(OUT_DIR, "greedy_sep_results.Rdata"))

# ── Load prior greedy paths ────────────────────────────────────────────────────

load(GREEDY_FILE)   # -> greedy_sep (list of 78 network results)

# Apply stopping rule to a stored path: return terms added before delta_aic >= 0
final_terms_from_path <- function(path_steps) {
  selected <- character(0)
  for (step in path_steps) {
    if (step$delta_aic >= 0) break
    selected <- c(selected, step$added)
  }
  selected
}

# Build lookup: file path -> (terms_init, terms_cont, cat, visible, synthetic)
sep_lookup <- setNames(lapply(greedy_sep, function(net) {
  list(
    terms_init = final_terms_from_path(net$path_init),
    terms_cont = final_terms_from_path(net$path_cont),
    cat        = net$cat,
    visible    = net$visible,
    synthetic  = if (!is.null(net$synthetic)) net$synthetic else NA
  )
}), sapply(greedy_sep, `[[`, "file"))

cat(sprintf("Loaded prior greedy paths for %d networks\n\n", length(sep_lookup)))

# ── Pooled candidate pool ──────────────────────────────────────────────────────
# Same pool as greedy_acc.R -- the only available baseline for pooled.

TERMS_POOL <- c(
  "SndInitiation", "RRecInitiation", "Repetition",
  "RecRecencyLog", "Imitation",      "IOviolation",
  "NOTPSnd",       "NITPSnd"
)

# Universe of terms that can appear in separable models (for coef extraction)
TERMS_SEP_INIT <- c(
  "RRecInitiation", "IOviolation", "InDegRec",  "NISPSnd",
  "NITPSnd",        "NOSPSnd",     "NOTPSnd",   "OutDegSnd"
)
TERMS_SEP_CONT <- c(
  "DyadRecencyLog", "IOviolation", "Imitation", "InDegRec",
  "NISPSnd",        "NITPSnd",     "NOSPSnd",   "NOTPSnd",
  "OutDegSnd",      "Repetition",  "WeightedRec","WeightedRep"
)

# ── Greedy forward selection (for pooled only) ─────────────────────────────────

greedy_forward <- function(pe, term_pool, min_events = 6L) {
  if (sum(pe$occurred) < min_events) return(NULL)
  null_fit <- tryCatch(
    coxph(Surv(time1, occurred) ~ strata(time2), data = pe),
    error = function(e) NULL
  )
  if (is.null(null_fit)) return(NULL)

  current_aic  <- AIC(null_fit)
  current_cols <- character(0)
  remaining    <- term_pool
  selected     <- character(0)

  repeat {
    if (length(remaining) == 0L) break
    candidates <- vapply(names(remaining), function(gname) {
      cols <- c(current_cols, remaining[[gname]])
      cols <- cols[cols %in% names(pe)]
      if (length(cols) == 0L) return(NA_real_)
      f <- as.formula(paste0(
        "Surv(time1, occurred) ~ ", paste(cols, collapse = " + "),
        " + strata(time2)"
      ))
      m <- tryCatch(coxph(f, data = pe), error = function(e) NULL)
      if (is.null(m)) return(NA_real_)
      a <- AIC(m)
      if (isTRUE(all.equal(a, current_aic, tolerance = 1e-6))) NA_real_ else a
    }, numeric(1))

    if (all(is.na(candidates))) break
    best_aic <- min(candidates, na.rm = TRUE)
    if (best_aic >= current_aic) break

    best_name    <- names(which.min(candidates))
    current_cols <- c(current_cols, remaining[[best_name]])
    selected     <- c(selected, best_name)
    remaining    <- remaining[names(remaining) != best_name]
    current_aic  <- best_aic
  }

  if (length(selected) == 0L) return(NULL)
  list(terms = selected, aic = current_aic, n_terms = length(selected))
}

# ── Predictive accuracy ────────────────────────────────────────────────────────

predict_accuracy <- function(fit, pe_scored, strata_col = "time2", regime = NULL) {
  if (is.null(fit)) return(NA_real_)
  cf <- tryCatch(coef(fit), error = function(e) NULL)
  if (is.null(cf) || all(is.na(cf))) return(NA_real_)
  cf <- cf[!is.na(cf)]
  tt <- tryCatch(model.matrix(fit, data = pe_scored), error = function(e) NULL)
  if (is.null(tt)) return(NA_real_)
  shared <- intersect(names(cf), colnames(tt))
  if (length(shared) == 0L) return(NA_real_)
  lp <- as.numeric(tt[, shared, drop = FALSE] %*% cf[shared])

  strata   <- pe_scored[[strata_col]]
  occurred <- pe_scored[["occurred"]]
  snd_init <- if (!is.null(regime)) pe_scored[["SndInitiation"]] else NULL

  correct <- vapply(unique(strata), function(s) {
    idx <- which(strata == s)
    obs <- which(occurred[idx] == 1L)
    if (length(obs) != 1L) return(NA_real_)
    if (!is.null(regime) && snd_init[idx[obs]] != regime) return(NA_real_)
    as.numeric(which.max(lp[idx]) == obs)
  }, numeric(1))

  mean(correct, na.rm = TRUE)
}

# ── Per-network function ───────────────────────────────────────────────────────

fit_one <- function(path, cat_label, visible_label, synthetic_label,
                    terms_init, terms_cont) {

  e <- new.env(parent = emptyenv())
  ok <- tryCatch({ load(path, envir = e); TRUE }, error = function(e) FALSE)
  if (!ok) return(NULL)

  el_raw <- tryCatch(
    as.data.frame(get("edgelist", envir = e))[, 1:3],
    error = function(e) NULL
  )
  if (is.null(el_raw)) return(NULL)
  el <- tryCatch(prepare_edgelist(el_raw), error = function(e) NULL)
  if (is.null(el) || nrow(el) < 10L) return(NULL)

  all_nodes <- sort(unique(c(el$sender, el$receiver)))
  base_rs   <- llrem:::get_base_riskset(NULL, all_nodes)
  pe <- tryCatch(
    llrem:::remstats_generate(el, n = length(all_nodes),
                              n_nonevents = 200L, base_riskset = base_rs),
    error = function(e) NULL
  )
  if (is.null(pe)) return(NULL)
  pe <- llrem:::add_survival_cols(pe)
  if (!"SndInitiation" %in% names(pe)) pe$SndInitiation <- pe$RecInitiation

  pe_init <- pe[pe$SndInitiation == 0L, ]
  pe_cont <- pe[pe$SndInitiation == 1L, ]
  n_init  <- sum(pe_init$occurred)
  n_cont  <- sum(pe_cont$occurred)
  if (n_init < 6L || n_cont < 6L) return(NULL)

  make_formula <- function(cols) {
    as.formula(paste0("Surv(time1, occurred) ~ ",
                      paste(cols, collapse = " + "),
                      " + strata(time2)"))
  }

  # ── Separable: refit the pre-determined final models ──────────────────────────

  fit_i <- if (length(terms_init) > 0L) {
    cols <- terms_init[terms_init %in% names(pe_init)]
    if (length(cols) > 0L)
      tryCatch(coxph(make_formula(cols), data = pe_init), error = function(e) NULL)
  }
  fit_c <- if (length(terms_cont) > 0L) {
    cols <- terms_cont[terms_cont %in% names(pe_cont)]
    if (length(cols) > 0L)
      tryCatch(coxph(make_formula(cols), data = pe_cont), error = function(e) NULL)
  }

  # ── Pooled: run greedy fresh ──────────────────────────────────────────────────

  pool_pool <- as.list(setNames(TERMS_POOL, TERMS_POOL))
  sel_pool  <- greedy_forward(pe, pool_pool)
  fit_p <- if (!is.null(sel_pool)) {
    cols <- sel_pool$terms[sel_pool$terms %in% names(pe)]
    tryCatch(coxph(make_formula(cols), data = pe), error = function(e) NULL)
  }

  # ── Coefficients ──────────────────────────────────────────────────────────────

  get_cf <- function(sm, nm) {
    if (!is.null(sm) && nm %in% rownames(sm)) sm[nm, "coef"]     else NA_real_
  }
  get_se <- function(sm, nm) {
    if (!is.null(sm) && nm %in% rownames(sm)) sm[nm, "se(coef)"] else NA_real_
  }
  get_z  <- function(sm, nm) {
    if (!is.null(sm) && nm %in% rownames(sm)) sm[nm, "z"]        else NA_real_
  }

  sm_i <- if (!is.null(fit_i)) tryCatch(coef(summary(fit_i)), error = function(e) NULL)
  sm_c <- if (!is.null(fit_c)) tryCatch(coef(summary(fit_c)), error = function(e) NULL)
  sm_p <- if (!is.null(fit_p)) tryCatch(coef(summary(fit_p)), error = function(e) NULL)

  # ── Predictive accuracy ───────────────────────────────────────────────────────

  acc_i <- predict_accuracy(fit_i, pe_init)
  acc_c <- predict_accuracy(fit_c, pe_cont)
  acc_sep <- if (!is.na(acc_i) && !is.na(acc_c))
    (acc_i * n_init + acc_c * n_cont) / (n_init + n_cont) else NA_real_

  acc_pool      <- predict_accuracy(fit_p, pe)
  acc_pool_init <- predict_accuracy(fit_p, pe, regime = 0L)
  acc_pool_cont <- predict_accuracy(fit_p, pe, regime = 1L)

  # ── Assemble ──────────────────────────────────────────────────────────────────

  out <- list(
    path          = path,
    cat           = cat_label,
    visible       = visible_label,
    synthetic     = synthetic_label,
    n             = length(all_nodes),
    n_events      = nrow(el),
    n_init        = n_init,
    n_cont        = n_cont,
    terms_init    = paste(terms_init, collapse = ","),
    terms_cont    = paste(terms_cont, collapse = ","),
    terms_pool    = if (!is.null(sel_pool)) paste(sel_pool$terms, collapse = ",") else NA_character_,
    n_terms_init  = length(terms_init),
    n_terms_cont  = length(terms_cont),
    n_terms_pool  = if (!is.null(sel_pool)) sel_pool$n_terms else NA_integer_,
    aic_sep       = if (!is.null(fit_i) && !is.null(fit_c))
                      fit_i$loglik[2] * -2 + 2 * length(coef(fit_i)) +
                      fit_c$loglik[2] * -2 + 2 * length(coef(fit_c)) else NA_real_,
    aic_pool      = if (!is.null(sel_pool)) sel_pool$aic else NA_real_,
    acc_sep_init  = acc_i,
    acc_sep_cont  = acc_c,
    acc_sep       = acc_sep,
    acc_pool      = acc_pool,
    acc_pool_init = acc_pool_init,
    acc_pool_cont = acc_pool_cont
  )

  for (nm in TERMS_SEP_INIT) {
    out[[paste0("init_", nm)]]        <- get_cf(sm_i, nm)
    out[[paste0("init_", nm, "_se")]] <- get_se(sm_i, nm)
    out[[paste0("init_", nm, "_z")]]  <- get_z(sm_i,  nm)
  }
  for (nm in TERMS_SEP_CONT) {
    out[[paste0("cont_", nm)]]        <- get_cf(sm_c, nm)
    out[[paste0("cont_", nm, "_se")]] <- get_se(sm_c, nm)
    out[[paste0("cont_", nm, "_z")]]  <- get_z(sm_c,  nm)
  }
  for (nm in TERMS_POOL) {
    out[[paste0("pool_", nm)]]        <- get_cf(sm_p, nm)
    out[[paste0("pool_", nm, "_se")]] <- get_se(sm_p, nm)
    out[[paste0("pool_", nm, "_z")]]  <- get_z(sm_p,  nm)
  }

  out
}

# ── Parallel execution ─────────────────────────────────────────────────────────

paths      <- names(sep_lookup)
cat_labels <- sapply(sep_lookup, `[[`, "cat")
vis_labels <- sapply(sep_lookup, `[[`, "visible")
syn_labels <- sapply(sep_lookup, `[[`, "synthetic")
ti_list    <- lapply(sep_lookup, `[[`, "terms_init")
tc_list    <- lapply(sep_lookup, `[[`, "terms_cont")

n_cores <- min(parallel::detectCores() - 1L, 6L)
cl <- parallel::makeCluster(n_cores)

parallel::clusterCall(cl, function() {
  devtools::load_all("llrem", quiet = TRUE)
  library(survival)
})

parallel::clusterExport(cl, c(
  "TERMS_POOL", "TERMS_SEP_INIT", "TERMS_SEP_CONT",
  "greedy_forward", "predict_accuracy", "fit_one", "DATA_DIR"
))

cat(sprintf("Running on %d networks using %d cores ...\n", length(paths), n_cores))

set.seed(8214L)
raw_results <- parallel::clusterMap(
  cl, fit_one,
  path             = paths,
  cat_label        = cat_labels,
  visible_label    = vis_labels,
  synthetic_label  = syn_labels,
  terms_init       = ti_list,
  terms_cont       = tc_list,
  SIMPLIFY         = FALSE
)
parallel::stopCluster(cl)

results <- Filter(Negate(is.null), raw_results)
cat(sprintf("%d networks retained (dropped %d)\n\n",
            length(results), length(paths) - length(results)))

res_df <- do.call(rbind, lapply(results, function(x) {
  as.data.frame(x[names(x) != "path"], stringsAsFactors = FALSE)
}))
rownames(res_df) <- sapply(results, `[[`, "path")

save(res_df, file = file.path(OUT_DIR, "greedy_acc2_results.Rdata"))
cat(sprintf("Saved to %s/greedy_acc2_results.Rdata\n\n", OUT_DIR))

# ── Summary ────────────────────────────────────────────────────────────────────

wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

n_per_domain <- table(res_df$cat)
res_df$w <- 1 / n_per_domain[res_df$cat]

cat("=== Mean predictive accuracy (unweighted) ===\n")
cat(sprintf("  Separable greedy (init):  %.4f\n", mean(res_df$acc_sep_init,  na.rm = TRUE)))
cat(sprintf("  Separable greedy (cont):  %.4f\n", mean(res_df$acc_sep_cont,  na.rm = TRUE)))
cat(sprintf("  Separable greedy (comb):  %.4f\n", mean(res_df$acc_sep,       na.rm = TRUE)))
cat(sprintf("  Pooled greedy:            %.4f\n", mean(res_df$acc_pool,      na.rm = TRUE)))
cat(sprintf("  Pooled greedy (init):     %.4f\n", mean(res_df$acc_pool_init, na.rm = TRUE)))
cat(sprintf("  Pooled greedy (cont):     %.4f\n", mean(res_df$acc_pool_cont, na.rm = TRUE)))

cat("\n=== Mean predictive accuracy (ecology-weighted) ===\n")
cat(sprintf("  Separable greedy (init):  %.4f\n", wmean(res_df$acc_sep_init,  res_df$w)))
cat(sprintf("  Separable greedy (cont):  %.4f\n", wmean(res_df$acc_sep_cont,  res_df$w)))
cat(sprintf("  Separable greedy (comb):  %.4f\n", wmean(res_df$acc_sep,       res_df$w)))
cat(sprintf("  Pooled greedy:            %.4f\n", wmean(res_df$acc_pool,      res_df$w)))
cat(sprintf("  Pooled greedy (init):     %.4f\n", wmean(res_df$acc_pool_init, res_df$w)))
cat(sprintf("  Pooled greedy (cont):     %.4f\n", wmean(res_df$acc_pool_cont, res_df$w)))

cat("\n=== Mean terms selected ===\n")
cat(sprintf("  Init sub-model:  %.1f\n", mean(res_df$n_terms_init, na.rm = TRUE)))
cat(sprintf("  Cont sub-model:  %.1f\n", mean(res_df$n_terms_cont, na.rm = TRUE)))
cat(sprintf("  Pooled:          %.1f / %d\n", mean(res_df$n_terms_pool, na.rm = TRUE), length(TERMS_POOL)))

cat("\nDone.\n")
