# covariate_demo.R
# Demonstrates the rem_covariates pipeline with simulated country-level data.
# Uses the international conflict edgelist bundled with llrem (1708 events,
# 95 nodes, time 1-1708).
#
# Covariates simulated:
#   polity      – static node:  regime type score (-10 autocracy … +10 democracy)
#   log_gdp     – temporal node: log GDP, two periods (pre/post t=850)
#   distance    – static dyad:  log geographic distance between country pairs
#   alliance    – temporal dyad: binary alliance membership, changes at t=500

library(llrem)

set.seed(2025)

# ── 1. Load and prepare edgelist ──────────────────────────────────────────────

data(edgelist)
el <- prepare_edgelist(edgelist)

n_nodes  <- length(unique(c(el$sender, el$receiver)))  # 95
all_nodes <- sort(unique(c(el$sender, el$receiver)))    # 1:95
n_events  <- nrow(el)                                   # 1708

cat(sprintf("Edgelist: %d events, %d nodes, time %d-%d\n",
            n_events, n_nodes, min(el$time), max(el$time)))

# ── 2. Simulate node-level covariates ─────────────────────────────────────────

# 2a. Static: Polity score (integer -10 to +10)
polity_df <- data.frame(
  node  = all_nodes,
  value = sample(-10:10, n_nodes, replace = TRUE)
)

# 2b. Temporal: log GDP, two time windows
#   Window 1: time  1–849  (early period)
#   Window 2: time 850–end (late period, modest growth ~+0.3 log units)
gdp_early <- rlnorm(n_nodes, meanlog = 10, sdlog = 1.5)
gdp_late  <- gdp_early * exp(rnorm(n_nodes, mean = 0.3, sd = 0.1))

log_gdp_df <- rbind(
  data.frame(node = all_nodes, time_start =    1, value = log(gdp_early)),
  data.frame(node = all_nodes, time_start =  850, value = log(gdp_late))
)

# ── 3. Simulate dyad-level covariates ─────────────────────────────────────────

all_dyads <- expand.grid(sender = all_nodes, receiver = all_nodes,
                          stringsAsFactors = FALSE)
all_dyads <- all_dyads[all_dyads$sender != all_dyads$receiver, ]

# 3a. Static: log geographic distance (symmetric, so simulate upper triangle
#     then mirror)
dist_mat <- matrix(NA_real_, n_nodes, n_nodes)
for (i in seq_len(n_nodes - 1L)) {
  for (j in (i + 1L):n_nodes) {
    d <- rlnorm(1, meanlog = 7, sdlog = 1.2)   # km, log-normal
    dist_mat[i, j] <- dist_mat[j, i] <- d
  }
}
distance_df <- data.frame(
  sender   = all_dyads$sender,
  receiver = all_dyads$receiver,
  value    = log(dist_mat[cbind(all_dyads$sender, all_dyads$receiver)])
)

# 3b. Temporal: binary alliance, two windows (t=1 and t=500 breakpoint)
#   ~15% of dyads allied in each period; alliance composition reshuffles

make_alliance_window <- function(dyads, prob_allied = 0.15) {
  # make it symmetric: if i-j allied then j-i allied
  n <- max(dyads$sender, dyads$receiver)
  allied_mat <- matrix(0L, n, n)
  pairs <- which(upper.tri(allied_mat), arr.ind = TRUE)
  allied_pairs <- pairs[sample(nrow(pairs), round(nrow(pairs) * prob_allied)), ]
  for (k in seq_len(nrow(allied_pairs))) {
    i <- allied_pairs[k, 1]; j <- allied_pairs[k, 2]
    allied_mat[i, j] <- allied_mat[j, i] <- 1L
  }
  data.frame(
    sender   = dyads$sender,
    receiver = dyads$receiver,
    value    = as.numeric(allied_mat[cbind(dyads$sender, dyads$receiver)])
  )
}

ally_early <- make_alliance_window(all_dyads, prob_allied = 0.12)
ally_late  <- make_alliance_window(all_dyads, prob_allied = 0.15)

alliance_df <- rbind(
  data.frame(ally_early, time_start =   1),
  data.frame(ally_late,  time_start = 500)
)

# ── 4. Build rem_covariates object ────────────────────────────────────────────

covs <- rem_covariates(
  # regime type: separate sender and receiver effects
  cov_node_static("polity", polity_df, role = "both"),
  # log GDP: separate sender and receiver effects, with temporal change
  cov_node_temporal("log_gdp", log_gdp_df, role = "both"),
  # log distance: dyad-level, static
  cov_dyad_static("log_dist", distance_df, warn_missing = FALSE),
  # alliance: dyad-level, temporal
  cov_dyad_temporal("alliance", alliance_df, warn_missing = FALSE)
)

print(covs)

# ── 5. Precompute the base risk set (cache for reuse) ─────────────────────────

base_rs <- llrem:::get_base_riskset(NULL, all_nodes)

# ── 6. Fit model with exogenous covariates ────────────────────────────────────
# Use first 500 events and clogit for speed.

el_sub <- el[seq_len(500L), ]

formula_with_covs <- survival::Surv(time1, occurred) ~
  Repetition + Recencylog + RecInitiation + IOviolation +
  NOTPSnd + NITPSnd + NOSPSnd + NISPSnd +
  polity_snd + polity_rec +
  log_gdp_snd + log_gdp_rec +
  log_dist + alliance +
  survival::strata(time2)

cat("\nFitting REM (clogit, 500 events, n_nonevents=200)...\n")
t_fit <- system.time(
  fit <- fit_rem(el_sub, formula_with_covs,
                 n            = n_nodes,
                 n_nonevents  = 200L,
                 base_riskset = base_rs,
                 method       = "clogit",
                 covariates   = covs)
)
cat(sprintf("Fit time: %.1fs\n\n", t_fit["elapsed"]))

# ── 7. Results ────────────────────────────────────────────────────────────────

print(fit)

# Quick check: any NAs in covariate columns of the risk set?
cat("\nChecking covariate NAs in risk set (n_nonevents=5 for speed)...\n")
pe <- llrem:::remstats_generate(el_sub,
                                 n            = n_nodes,
                                 n_nonevents  = 5L,
                                 base_riskset = base_rs,
                                 covariates   = covs)
cov_cols <- c("polity_snd", "polity_rec", "log_gdp_snd", "log_gdp_rec",
               "log_dist", "alliance")
na_counts <- colSums(is.na(pe[, cov_cols]))
cat("NA counts per covariate column:\n")
print(na_counts)
cat("\nEndogenous NA check (should all be 0):\n")
print(colSums(is.na(pe[, c("Repetition", "Recencylog", "NTriadic")])))
