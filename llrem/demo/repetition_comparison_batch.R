# repetition_comparison_batch.R
# Fits Repetition vs RSndSnd models (llrem clogit only) across all DataP
# datasets and reports delta log-likelihood for each.
#
# Uses llrem only — relevent is too slow at full risk set to loop across
# datasets interactively.  The llrem result is sufficient to test whether
# RSndSnd consistently out-fits Repetition.

library(llrem)

set.seed(2025)

N_NONEVENTS <- 200L
N_MAX       <- 1000L  # cap events per dataset

DATA_DIR <- Sys.getenv("DATA_DIR", "")   # path to the DataP network data directory

datasets <- list(
  Classroom = file.path(DATA_DIR, "Classroom/McFarland.Rdata"),
  TSOF      = file.path(DATA_DIR, "Email/TSOF.Rdata"),
  EU1       = file.path(DATA_DIR, "Email/eu1.Rdata"),
  EU2       = file.path(DATA_DIR, "Email/eu2.Rdata"),
  EU3       = file.path(DATA_DIR, "Email/eu3.Rdata"),
  EU4       = file.path(DATA_DIR, "Email/eu4.Rdata"),
  UCIMsg    = file.path(DATA_DIR, "PMs/ucimsg.Rdata"),
  Phone     = file.path(DATA_DIR, "Phone/Phone.Rdata"),
  CrimePhone= file.path(DATA_DIR, "Phone/crimephone.Rdata"),
  SMS       = file.path(DATA_DIR, "SMS/sms1.Rdata"),
  Twitter   = file.path(DATA_DIR, "Twitter/twitter_boston_relevent.rdata"),
  MID       = "bundled"
)

fit_pair <- function(el_sub, n_nodes, base_rs) {
  f_rep  <- survival::Surv(time1, occurred) ~
    Repetition + RecInitiation + Recencylog + IOviolation + survival::strata(time2)
  f_rsnd <- survival::Surv(time1, occurred) ~
    RSndSnd    + RecInitiation + Recencylog + IOviolation + survival::strata(time2)

  fa <- tryCatch(
    fit_rem(el_sub, f_rep,  n=n_nodes, n_nonevents=N_NONEVENTS,
            base_riskset=base_rs, method="clogit"),
    error = function(e) NULL)
  fb <- tryCatch(
    fit_rem(el_sub, f_rsnd, n=n_nodes, n_nonevents=N_NONEVENTS,
            base_riskset=base_rs, method="clogit"),
    error = function(e) NULL)
  list(rep=fa, rsnd=fb)
}

results <- list()

for (nm in names(datasets)) {
  path <- datasets[[nm]]
  cat(sprintf("%-12s  ", nm))

  if (path == "bundled") {
    data(edgelist, package="llrem")
    el_raw <- edgelist
    n_nodes <- length(unique(c(el_raw$id_snd, el_raw$id_rec)))
  } else {
    e <- new.env(parent=emptyenv())
    load(path, envir=e)
    el_raw  <- as.data.frame(get("edgelist", envir=e))
    n_nodes <- if (exists("n", envir=e)) as.integer(get("n", envir=e))
               else length(unique(c(el_raw[[2]], el_raw[[3]])))
  }

  el      <- prepare_edgelist(el_raw)
  n_ev    <- min(nrow(el), N_MAX)
  el_sub  <- el[seq_len(n_ev), ]
  all_nodes <- sort(unique(c(el$sender, el$receiver)))
  # use full node range so IDs outside the subset are valid
  base_rs <- llrem:::get_base_riskset(NULL, seq_len(max(all_nodes)))

  fits <- fit_pair(el_sub, max(all_nodes), base_rs)

  if (is.null(fits$rep) || is.null(fits$rsnd)) {
    cat("ERROR\n")
    next
  }

  ll_rep  <- fits$rep$loglik   # c(null_negll, fitted_negll)
  ll_rsnd <- fits$rsnd$loglik
  delta   <- ll_rsnd[2] - ll_rep[2]   # positive = RSndSnd fits better

  results[[nm]] <- list(
    n_nodes = max(all_nodes),
    n_events = n_ev,
    ll_rep   = ll_rep[2],
    ll_rsnd  = ll_rsnd[2],
    delta    = delta,
    converged_rep  = fits$rep$converged,
    converged_rsnd = fits$rsnd$converged
  )

  cat(sprintf("n=%3d  m=%4d  llRep=%9.2f  llRSnd=%9.2f  delta=%+7.2f  %s\n",
              max(all_nodes), n_ev,
              ll_rep[2], ll_rsnd[2], delta,
              if (delta > 0) "RSndSnd better" else "Repetition better"))
}

# ── Summary ───────────────────────────────────────────────────────────────────

cat("\n")
cat(strrep("=", 65), "\n")
cat("Summary: delta loglik = ll(RSndSnd) - ll(Repetition)\n")
cat("Positive = RSndSnd fits better (lower negll = better clogit fit)\n")
cat(strrep("=", 65), "\n")

deltas <- sapply(results, `[[`, "delta")
cat(sprintf("\n  RSndSnd better in %d / %d datasets\n",
            sum(deltas > 0), length(deltas)))
cat(sprintf("  Median delta: %+.2f\n", median(deltas)))
cat(sprintf("  Range:        %+.2f to %+.2f\n", min(deltas), max(deltas)))
