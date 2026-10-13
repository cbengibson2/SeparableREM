# relevent_comparison.R
# Coefficient comparison: llrem (clogit) vs relevent::rem.dyad(ordinal=TRUE)
# on the first 200 events of the bundled MID conflict edgelist.
#
# Settled model specification:
#   llrem             relevent          Description
#   -----------       --------          -----------
#   Repetition        FrPSndSnd         dyad_occ / sender_outDeg
#   RecInitiation     (none)            has i→j ever fired? (dyad persistence)
#   Recencylog        RRecSnd           reverse-dyad recency (reciprocity)
#   RRecInitiation    PSAB-BA           has j→i ever fired? (reciprocity binary)
#   IOviolation       PSAB-XA/XB/XY     interaction-order p-shifts (collapsed vs separate)
#
# Note: RRecInitiation captures the binary threshold of Recencylog (ever vs never);
# in relevent's p-shift framework the closest analogue is PSAB-BA (turn-response).

library(llrem)
library(relevent)

set.seed(2025)

# ── 1. Data ───────────────────────────────────────────────────────────────────

data(edgelist)
el <- prepare_edgelist(edgelist)

# 200 events: relevent's ordinal likelihood is O(m * n^2) with no sampling;
# at n=95, 200 events ≈ 90s.
N_EVENTS  <- 500L
el_sub    <- el[seq_len(N_EVENTS), ]
el_mat    <- as.matrix(el_sub[, c("time", "sender", "receiver")])

n_nodes   <- length(unique(c(el$sender, el$receiver)))
all_nodes <- sort(unique(c(el$sender, el$receiver)))
base_rs   <- llrem:::get_base_riskset(NULL, all_nodes)

# ── 2. llrem clogit ───────────────────────────────────────────────────────────

cat(sprintf("Fitting llrem clogit (%d events, n_nonevents=200)...\n", N_EVENTS))
t1 <- system.time(
  fit_llrem <- fit_rem(el_sub,
    survival::Surv(time1, occurred) ~
      Repetition + RecInitiation + Recencylog + RRecInitiation + IOviolation +
      survival::strata(time2),
    n = n_nodes, n_nonevents = 200L,
    base_riskset = base_rs, method = "clogit")
)
cat(sprintf("  Done (%.1fs)\n\n", t1["elapsed"]))

# ── 3. relevent rem.dyad ordinal ──────────────────────────────────────────────

cat(sprintf("Fitting relevent rem.dyad (ordinal=TRUE, MLE, %d events)...\n", N_EVENTS))
cat("  (full risk set; ~90s at n=95)\n")
t2 <- system.time(
  fit_rel <- rem.dyad(el_mat, n = n_nodes,
    effects    = c("FrPSndSnd", "RRecSnd", "PSAB-BA",
                   "PSAB-XA", "PSAB-XB", "PSAB-XY"),
    ordinal    = TRUE, fit.method = "MLE",
    hessian    = TRUE, verbose = FALSE, gof = FALSE)
)
cat(sprintf("  Done (%.1fs)\n\n", t2["elapsed"]))

# ── 4. Coefficient table ──────────────────────────────────────────────────────

ll_coef <- fit_llrem$coef
ll_se   <- sqrt(diag(fit_llrem$vcov))

rel_coef <- fit_rel$coef
rel_se   <- tryCatch(
  if (!is.null(fit_rel$cov)) sqrt(diag(fit_rel$cov))
  else sqrt(diag(solve(fit_rel$hessian / 2))),
  error = function(e) rep(NA_real_, length(fit_rel$coef))
)
names(rel_se) <- names(rel_coef)

cat("=============================================================\n")
cat("Coefficient comparison: llrem (clogit) vs relevent (ordinal)\n")
cat(sprintf("First %d MID events | llrem n_nonevents=200 | relevent full risk set\n",
            N_EVENTS))
cat("=============================================================\n\n")

cat(sprintf("%-18s  %8s %8s  |  %-20s %8s %8s\n",
            "llrem", "coef", "se", "relevent", "coef", "se"))
cat(strrep("-", 75), "\n")

# Simple paired rows
simple_pairs <- list(
  c("Repetition",     "FrPSndSnd"),
  c("Recencylog",     "RRecSnd"),
  c("RRecInitiation", "PSAB-BA")
)
for (pr in simple_pairs)
  cat(sprintf("%-18s  %8.4f %8.4f  |  %-20s %8.4f %8.4f\n",
              pr[1], ll_coef[pr[1]], ll_se[pr[1]],
              pr[2], rel_coef[pr[2]], rel_se[pr[2]]))

# RecInitiation: no relevent equivalent
cat(sprintf("%-18s  %8.4f %8.4f  |  %-20s\n",
            "RecInitiation", ll_coef["RecInitiation"], ll_se["RecInitiation"],
            "(no equivalent)"))

# IOviolation: three separate relevent p-shifts
pshifts <- c("PSAB-XA", "PSAB-XB", "PSAB-XY")
cat(sprintf("%-18s  %8.4f %8.4f  |  %-20s %8.4f %8.4f\n",
            "IOviolation",
            ll_coef["IOviolation"], ll_se["IOviolation"],
            pshifts[1], rel_coef[pshifts[1]], rel_se[pshifts[1]]))
for (ps in pshifts[-1])
  cat(sprintf("%-18s  %8s %8s  |  %-20s %8.4f %8.4f\n",
              "", "", "", ps, rel_coef[ps], rel_se[ps]))

cat(strrep("-", 75), "\n\n")

# ── 5. Log-likelihood comparison ──────────────────────────────────────────────

rel_null_ll <- -(N_EVENTS - 1) * log(n_nodes * (n_nodes - 1L))
cat(sprintf("llrem    loglik: %.2f  (null: %.2f)  dLL: %.2f\n",
            fit_llrem$loglik[2], fit_llrem$loglik[1],
            fit_llrem$loglik[2] - fit_llrem$loglik[1]))
cat(sprintf("relevent loglik: %.2f  (null: %.2f)  dLL: %.2f\n\n",
            -fit_rel$value, rel_null_ll, -fit_rel$value - rel_null_ll))

# ── 6. Notes ──────────────────────────────────────────────────────────────────

cat("Notes:\n")
cat("  RecInitiation  : binary — has i→j ever fired? (dyad persistence)\n")
cat("    No relevent equivalent; FrPSndSnd absorbs this unconditionally.\n\n")
cat("  RRecInitiation : binary — has j→i ever fired? (reciprocity initiation)\n")
cat("    Captures threshold component of Recencylog (ever vs never).\n")
cat("    relevent PSAB-BA is the sequential analogue: B responds to A→B with B→A.\n\n")
cat("  Recencylog sign: −log1p(t) so negative coef = more-recent reverse tie\n")
cat("    raises rate. RRecSnd is a rank so positive coef = same direction.\n\n")
cat("  IOviolation collapses PSAB-XA + PSAB-XB + PSAB-XY into one binary.\n")
