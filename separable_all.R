# separable_all.R
#
# Separable relational event models across all available datasets.
#
# Estimation: two separate Cox models on disjoint event subsets.
#
# Initiation model  (SndInitiation == 0, i has never previously sent to j):
#   RRecInitiation, NOTPSnd, NITPSnd, Recencylog, Imitation, IOviolation
#
# Continuation model (SndInitiation == 1, i has sent to j before):
#   Repetition, NOTPSnd, NITPSnd, Recencylog, Imitation, IOviolation
#
# Pooled reference model (all events, shared coefficients):
#   SndInitiation + RRecInitiation + Repetition + NOTPSnd + NITPSnd +
#   Recencylog + Imitation + IOviolation
#
# Separable df = 6 + 6 = 12; pooled df = 8.
# Model comparison uses AIC: delta_AIC = AIC(pool) - AIC(sep),
# positive = separable preferred.

rm(list = ls())
devtools::load_all("llrem", quiet = TRUE)
library(parallel)
library(survival)

# ── Domain metadata ───────────────────────────────────────────────────────────

DATA_DIR <- Sys.getenv("DATA_DIR", "")   # path to the DataP network data directory

domain_meta <- list(
  Classroom = list(dir = "Classroom", cat = "Classroom",    visible = "Y"),
  Email     = list(dir = "Email",     cat = "Email",        visible = "N"),
  Gangs     = list(dir = "Gangs",     cat = "Gang Violence",visible = "Y"),
  PMs       = list(dir = "PMs",       cat = "PMs",          visible = "N"),
  Phone     = list(dir = "Phone",     cat = "Phone (DoD)",  visible = "Y"),
  SMS       = list(dir = "SMS",       cat = "SMS",          visible = "N"),
  Twitter   = list(dir = "Twitter",   cat = "Twitter",      visible = "Y"),
  War       = list(dir = "War",       cat = "War",          visible = "Y"),
  ftf       = list(dir = "ftf",       cat = "Glances",      visible = "N")
)

all_files <- do.call(rbind, lapply(names(domain_meta), function(dname) {
  meta  <- domain_meta[[dname]]
  flist <- list.files(file.path(DATA_DIR, meta$dir),
                      pattern = "\\.[Rr]data$", full.names = TRUE)
  data.frame(path = flist, cat = meta$cat, visible = meta$visible,
             stringsAsFactors = FALSE)
}))

cat(sprintf("Found %d files across %d domains\n\n",
            nrow(all_files), length(unique(all_files$cat))))
print(table(all_files$cat))

# ── Formulas ──────────────────────────────────────────────────────────────────

form_init <- Surv(time1, occurred) ~
  RRecInitiation + NOTPSnd + NITPSnd + Recencylog + Imitation + IOviolation +
  strata(time2)

form_cont <- Surv(time1, occurred) ~
  Repetition + NOTPSnd + NITPSnd + Recencylog + Imitation + IOviolation +
  strata(time2)

form_pool <- Surv(time1, occurred) ~
  SndInitiation + RRecInitiation + Repetition + NOTPSnd + NITPSnd +
  Recencylog + Imitation + IOviolation +
  strata(time2)

# ── Per-network fitting function ──────────────────────────────────────────────

fit_one <- function(path, cat_label, visible_label) {
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
  pe$SndInitiation <- pe$RecInitiation

  pe_init <- pe[pe$SndInitiation == 0, ]
  pe_cont <- pe[pe$SndInitiation == 1, ]

  n_init_ev <- sum(pe_init$occurred)
  n_cont_ev <- sum(pe_cont$occurred)

  # Need at least as many events as parameters in each sub-model
  if (n_init_ev < 6L || n_cont_ev < 6L) return(NULL)

  fit_i <- tryCatch(coxph(form_init, data = pe_init), error = function(e) NULL)
  fit_c <- tryCatch(coxph(form_cont, data = pe_cont), error = function(e) NULL)
  fit_p <- tryCatch(coxph(form_pool, data = pe),      error = function(e) NULL)

  if (is.null(fit_i) || is.null(fit_c) || is.null(fit_p)) return(NULL)

  cf_i <- coef(summary(fit_i))
  cf_c <- coef(summary(fit_c))
  cf_p <- coef(summary(fit_p))

  get_cf <- function(cf, nm, col = "coef") {
    if (nm %in% rownames(cf)) cf[nm, col] else NA_real_
  }
  get_z <- function(cf, nm) get_cf(cf, nm, "z")

  ll_i   <- fit_i$loglik[2]
  ll_c   <- fit_c$loglik[2]
  ll_p   <- fit_p$loglik[2]
  ll_sep <- ll_i + ll_c

  df_i   <- length(coef(fit_i))   # 6
  df_c   <- length(coef(fit_c))   # 6
  df_p   <- length(coef(fit_p))   # 8
  df_sep <- df_i + df_c

  aic_sep  <- -2 * ll_sep + 2 * df_sep
  aic_pool <- -2 * ll_p   + 2 * df_p
  delta_aic <- aic_pool - aic_sep   # positive = separable preferred

  list(
    path     = path,
    cat      = cat_label,
    visible  = visible_label,
    n        = length(all_nodes),
    n_events = nrow(el),
    n_init   = n_init_ev,
    n_cont   = n_cont_ev,

    # ── Initiation model coefficients ─────────────────────────────────────────
    init_RRecInitiation   = get_cf(cf_i, "RRecInitiation"),
    init_NOTPSnd          = get_cf(cf_i, "NOTPSnd"),
    init_NITPSnd          = get_cf(cf_i, "NITPSnd"),
    init_Recencylog       = get_cf(cf_i, "Recencylog"),
    init_Imitation        = get_cf(cf_i, "Imitation"),
    init_IOviolation      = get_cf(cf_i, "IOviolation"),

    init_RRecInitiation_z = get_z(cf_i, "RRecInitiation"),
    init_NOTPSnd_z        = get_z(cf_i, "NOTPSnd"),
    init_NITPSnd_z        = get_z(cf_i, "NITPSnd"),
    init_Recencylog_z     = get_z(cf_i, "Recencylog"),
    init_Imitation_z      = get_z(cf_i, "Imitation"),
    init_IOviolation_z    = get_z(cf_i, "IOviolation"),

    # ── Continuation model coefficients ───────────────────────────────────────
    cont_Repetition       = get_cf(cf_c, "Repetition"),
    cont_NOTPSnd          = get_cf(cf_c, "NOTPSnd"),
    cont_NITPSnd          = get_cf(cf_c, "NITPSnd"),
    cont_Recencylog       = get_cf(cf_c, "Recencylog"),
    cont_Imitation        = get_cf(cf_c, "Imitation"),
    cont_IOviolation      = get_cf(cf_c, "IOviolation"),

    cont_Repetition_z     = get_z(cf_c, "Repetition"),
    cont_NOTPSnd_z        = get_z(cf_c, "NOTPSnd"),
    cont_NITPSnd_z        = get_z(cf_c, "NITPSnd"),
    cont_Recencylog_z     = get_z(cf_c, "Recencylog"),
    cont_Imitation_z      = get_z(cf_c, "Imitation"),
    cont_IOviolation_z    = get_z(cf_c, "IOviolation"),

    # ── Pooled model coefficients ──────────────────────────────────────────────
    pool_SndInitiation    = get_cf(cf_p, "SndInitiation"),
    pool_RRecInitiation   = get_cf(cf_p, "RRecInitiation"),
    pool_Repetition       = get_cf(cf_p, "Repetition"),
    pool_NOTPSnd          = get_cf(cf_p, "NOTPSnd"),
    pool_NITPSnd          = get_cf(cf_p, "NITPSnd"),
    pool_Recencylog       = get_cf(cf_p, "Recencylog"),
    pool_Imitation        = get_cf(cf_p, "Imitation"),
    pool_IOviolation      = get_cf(cf_p, "IOviolation"),

    # ── Model fit ──────────────────────────────────────────────────────────────
    ll_init   = ll_i,
    ll_cont   = ll_c,
    ll_sep    = ll_sep,
    ll_pool   = ll_p,
    df_sep    = df_sep,
    df_pool   = df_p,
    aic_sep   = aic_sep,
    aic_pool  = aic_pool,
    delta_aic = delta_aic
  )
}

# ── Parallel execution ────────────────────────────────────────────────────────

n_cores <- min(parallel::detectCores() - 1L, 4L)
cl <- parallel::makeCluster(n_cores)
parallel::clusterCall(cl, function() {
  devtools::load_all("llrem", quiet = TRUE)
  library(survival)
})
parallel::clusterExport(cl, c("form_init", "form_cont", "form_pool",
                               "fit_one", "DATA_DIR"))

cat(sprintf("\nFitting separable models on %d networks using %d cores ...\n",
            nrow(all_files), n_cores))

set.seed(8214L)
raw_results <- parallel::clusterMap(
  cl, fit_one,
  path          = all_files$path,
  cat_label     = all_files$cat,
  visible_label = all_files$visible,
  SIMPLIFY      = FALSE
)
parallel::stopCluster(cl)

results <- raw_results[!sapply(raw_results, is.null)]
cat(sprintf("%d networks retained (dropped %d)\n\n",
            length(results), nrow(all_files) - length(results)))

# ── Flatten to data frame ─────────────────────────────────────────────────────

res_df <- do.call(rbind, lapply(results, function(x) {
  as.data.frame(x[names(x) != "path"], stringsAsFactors = FALSE)
}))
res_df$vis_label <- ifelse(res_df$visible == "Y", "High", "Low")

save(res_df, file = "separable_all_results.Rdata")
cat("Raw results saved to separable_all_results.Rdata\n\n")

# ── Summary helpers ───────────────────────────────────────────────────────────

smry <- function(x, label = "") {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(data.frame(label=label, n=0, mean=NA, sd=NA, t=NA))
  m <- mean(x); s <- sd(x); n <- length(x)
  data.frame(label = label, n = n,
             mean  = round(m, 3),
             sd    = round(s, 3),
             t     = round(m / (s / sqrt(n)), 2))
}

smry_by_cat <- function(col, df, label) {
  do.call(rbind, lapply(sort(unique(df$cat)), function(g) {
    r <- smry(df[[col]][df$cat == g], label = g)
    r$term <- label
    r
  }))
}

# ── Table 1: Initiation model coefficients by domain ─────────────────────────

cat("=================================================================\n")
cat("TABLE 1: INITIATION MODEL — mean coefficient by domain\n")
cat("=================================================================\n")
init_terms <- c("init_RRecInitiation","init_NOTPSnd","init_NITPSnd",
                "init_Recencylog","init_Imitation","init_IOviolation")
init_labs  <- c("RRecInitiation","NOTPSnd","NITPSnd",
                "Recencylog","Imitation","IOviolation")

for (i in seq_along(init_terms)) {
  cat(sprintf("\n  --- %s ---\n", init_labs[i]))
  print(smry_by_cat(init_terms[i], res_df, init_labs[i])[, c("label","n","mean","sd","t")],
        row.names = FALSE)
}

# ── Table 2: Continuation model coefficients by domain ───────────────────────

cat("\n=================================================================\n")
cat("TABLE 2: CONTINUATION MODEL — mean coefficient by domain\n")
cat("=================================================================\n")
cont_terms <- c("cont_Repetition","cont_NOTPSnd","cont_NITPSnd",
                "cont_Recencylog","cont_Imitation","cont_IOviolation")
cont_labs  <- c("Repetition","NOTPSnd","NITPSnd",
                "Recencylog","Imitation","IOviolation")

for (i in seq_along(cont_terms)) {
  cat(sprintf("\n  --- %s ---\n", cont_labs[i]))
  print(smry_by_cat(cont_terms[i], res_df, cont_labs[i])[, c("label","n","mean","sd","t")],
        row.names = FALSE)
}

# ── Table 3: Pooled model coefficients by domain ─────────────────────────────

cat("\n=================================================================\n")
cat("TABLE 3: POOLED MODEL — mean coefficient by domain\n")
cat("         (compare sign/magnitude to Tables 1 and 2)\n")
cat("=================================================================\n")
pool_terms <- c("pool_SndInitiation","pool_RRecInitiation","pool_Repetition",
                "pool_NOTPSnd","pool_NITPSnd","pool_Recencylog",
                "pool_Imitation","pool_IOviolation")
pool_labs  <- c("SndInitiation","RRecInitiation","Repetition",
                "NOTPSnd","NITPSnd","Recencylog","Imitation","IOviolation")

for (i in seq_along(pool_terms)) {
  cat(sprintf("\n  --- %s ---\n", pool_labs[i]))
  print(smry_by_cat(pool_terms[i], res_df, pool_labs[i])[, c("label","n","mean","sd","t")],
        row.names = FALSE)
}

# ── Table 4: AIC comparison by domain ────────────────────────────────────────

cat("\n=================================================================\n")
cat("TABLE 4: MODEL FIT — separable vs pooled, by domain\n")
cat("  delta_AIC = AIC(pool) - AIC(sep)  (positive = separable preferred)\n")
cat(sprintf("  df_sep = %d (init 6 + cont 6),  df_pool = 8\n",
            12L))
cat("=================================================================\n")

fit_tab <- do.call(rbind, lapply(sort(unique(res_df$cat)), function(g) {
  sub <- res_df[res_df$cat == g, ]
  data.frame(
    Domain    = g,
    N         = nrow(sub),
    ll_sep    = round(mean(sub$ll_sep,    na.rm = TRUE), 1),
    ll_pool   = round(mean(sub$ll_pool,   na.rm = TRUE), 1),
    aic_sep   = round(mean(sub$aic_sep,   na.rm = TRUE), 1),
    aic_pool  = round(mean(sub$aic_pool,  na.rm = TRUE), 1),
    delta_aic = round(mean(sub$delta_aic, na.rm = TRUE), 1),
    pct_pos   = round(100 * mean(sub$delta_aic > 0, na.rm = TRUE), 1)
  )
}))
print(fit_tab, row.names = FALSE)

cat(sprintf("\nOverall: separable preferred (delta_AIC > 0) in %.1f%% of networks\n",
            100 * mean(res_df$delta_aic > 0, na.rm = TRUE)))

# ── Table 5: RRecInitiation — initiation vs pooled (sign flip diagnostic) ────

cat("\n=================================================================\n")
cat("TABLE 5: RRecInitiation — initiation model vs pooled\n")
cat("  Sign flip indicates pooled model confounds initiation/continuation\n")
cat("=================================================================\n")

flip_tab <- do.call(rbind, lapply(sort(unique(res_df$cat)), function(g) {
  sub <- res_df[res_df$cat == g, ]
  data.frame(
    Domain        = g,
    N             = nrow(sub),
    init_mean     = round(mean(sub$init_RRecInitiation,  na.rm = TRUE), 3),
    pool_mean     = round(mean(sub$pool_RRecInitiation,  na.rm = TRUE), 3),
    sign_flip_pct = round(100 * mean(
      sign(sub$init_RRecInitiation) != sign(sub$pool_RRecInitiation),
      na.rm = TRUE), 1)
  )
}))
print(flip_tab, row.names = FALSE)

# ── Table 6: NOTPSnd — initiation vs continuation vs pooled ──────────────────

cat("\n=================================================================\n")
cat("TABLE 6: NOTPSnd — initiation vs continuation vs pooled\n")
cat("  Tests whether transitive closure is initiation-specific\n")
cat("=================================================================\n")

notp_tab <- do.call(rbind, lapply(sort(unique(res_df$cat)), function(g) {
  sub <- res_df[res_df$cat == g, ]
  data.frame(
    Domain    = g,
    N         = nrow(sub),
    init_mean = round(mean(sub$init_NOTPSnd,  na.rm = TRUE), 3),
    cont_mean = round(mean(sub$cont_NOTPSnd,  na.rm = TRUE), 3),
    pool_mean = round(mean(sub$pool_NOTPSnd,  na.rm = TRUE), 3)
  )
}))
print(notp_tab, row.names = FALSE)

# ── Table 7: Split by visibility ─────────────────────────────────────────────

cat("\n=================================================================\n")
cat("TABLE 7: KEY COEFFICIENTS BY VISIBILITY GROUP\n")
cat("=================================================================\n")

vis_tab <- do.call(rbind, lapply(c("High","Low"), function(v) {
  sub <- res_df[res_df$vis_label == v, ]
  data.frame(
    Visibility         = v,
    N                  = nrow(sub),
    init_RRecInitiation = round(mean(sub$init_RRecInitiation, na.rm=TRUE), 3),
    init_NOTPSnd        = round(mean(sub$init_NOTPSnd,        na.rm=TRUE), 3),
    init_NITPSnd        = round(mean(sub$init_NITPSnd,        na.rm=TRUE), 3),
    cont_Repetition     = round(mean(sub$cont_Repetition,     na.rm=TRUE), 3),
    cont_Recencylog     = round(mean(sub$cont_Recencylog,     na.rm=TRUE), 3),
    delta_aic           = round(mean(sub$delta_aic,           na.rm=TRUE), 2)
  )
}))
print(vis_tab, row.names = FALSE)

# ── Table 8: Event partition — initiation fraction by domain ─────────────────

cat("\n=================================================================\n")
cat("TABLE 8: EVENT PARTITION — initiation fraction by domain\n")
cat("=================================================================\n")

part_tab <- do.call(rbind, lapply(sort(unique(res_df$cat)), function(g) {
  sub <- res_df[res_df$cat == g, ]
  data.frame(
    Domain      = g,
    N_networks  = nrow(sub),
    mean_events = round(mean(sub$n_events, na.rm=TRUE)),
    mean_init   = round(mean(sub$n_init,   na.rm=TRUE)),
    mean_cont   = round(mean(sub$n_cont,   na.rm=TRUE)),
    init_pct    = round(100 * mean(sub$n_init / sub$n_events, na.rm=TRUE), 1)
  )
}))
print(part_tab, row.names = FALSE)

cat("\nDone.\n")
