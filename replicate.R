# replicate.R
#
# Replication script for all numbers in the paper.
# Reads from authoritative stored results; does not re-fit any models.
# Prints every claim from the paper (section by section) with the
# computed value alongside it so deviations are immediately visible.
#
# Run order:
#   1. revision.R            -> revision_results.Rdata   (primary models)
#   2. greedy_opt_sep.R      -> greedy_sep_results.Rdata (greedy paths)
#   3. greedy_acc2.R         -> greedy_acc2_results.Rdata (greedy accuracy)
#   4. This script           -> console audit only
#
# All hard-coded comparators are taken directly from the paper text / tables.

rm(list = ls())

load("revision_results.Rdata")
load("greedy_sep_results.Rdata")
load("greedy_acc2_results.Rdata")

domain_order <- c("Glances", "Email", "War", "SMS", "Classroom",
                  "Twitter", "Phone (DoD)", "PMs", "Gang Violence")

# ── Helpers ───────────────────────────────────────────────────────────────────

agg <- function(sub, bcol, zcol) {
  b <- sub[[bcol]][is.finite(sub[[bcol]])]
  z <- sub[[zcol]][is.finite(sub[[zcol]])]
  N <- length(b)
  if (N == 0) return(list(mean_b = NA_real_, tstat = NA_real_, N = 0L))
  m <- mean(b)
  t <- if (N == 1) { if (length(z) > 0) z[1] else NA_real_ } else {
    s <- sd(b); if (!is.na(s) && s > 0) m / (s / sqrt(N)) else NA_real_
  }
  list(mean_b = m, tstat = t, N = N)
}

domain_agg <- function(bcol, zcol) {
  lapply(setNames(domain_order, domain_order), function(d) {
    agg(res_df[res_df$cat == d, ], bcol, zcol)
  })
}

chk <- function(label, computed, paper, tol = 0.005) {
  ok  <- !is.na(computed) && !is.na(paper) && abs(computed - paper) <= tol
  flag <- if (is.na(computed)) "  [MISSING]" else
          if (ok)              "  OK" else
          sprintf("  MISMATCH: paper=%.4f  computed=%.4f  diff=%.4f",
                  paper, computed, computed - paper)
  cat(sprintf("  %-65s%s\n", label, flag))
}

section <- function(s) {
  cat("\n", strrep("=", 72), "\n", s, "\n", strrep("=", 72), "\n", sep = "")
}

# ── Table 1: Event partition (tab:partition) ──────────────────────────────────
section("Table 1: Event partition (tab:partition)")

part <- do.call(rbind, lapply(domain_order, function(d) {
  sub <- res_df[res_df$cat == d, ]
  data.frame(
    domain    = d,
    N         = nrow(sub),
    mean_ev   = round(mean(sub$n_events)),
    mean_init = round(mean(sub$n_init)),
    init_pct  = round(100 * mean(sub$n_init / sub$n_events), 1),
    stringsAsFactors = FALSE
  )
}))

paper_part <- data.frame(
  domain    = domain_order,
  N         = c(62, 33,   1,   3,   1,   1,   1,   8, 100),
  mean_ev   = c(10000, 5145, 1708, 5001, 691, 227, 196, 7001, 450),
  mean_init = c(43, 446, 137, 589, 113,  48,  51, 2966, 209),
  init_pct  = c(0.4, 8.6, 8.0, 11.8, 16.4, 21.1, 26.0, 42.4, 46.4),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(part))) {
  d <- part$domain[i]
  chk(paste0(d, " N"),         part$N[i],         paper_part$N[i],         tol = 0)
  chk(paste0(d, " init_pct"),  part$init_pct[i],  paper_part$init_pct[i],  tol = 0.2)
}

# ── Text: initiation fraction section ────────────────────────────────────────
section("Text: initiation fraction section")

gv_sub  <- res_df[res_df$cat == "Gang Violence", ]
gl_sub  <- res_df[res_df$cat == "Glances", ]
rep_gv  <- agg(gv_sub, "cont_Repetition", "cont_Repetition_z")

chk("Gang Violence actors = 33",
    mean(sqrt(gv_sub$n)),     # n is n_actors^2-ish -- check n field directly
    NA)  # not directly computable from res_df; skip
chk("Gang Violence mean events = 450",  mean(gv_sub$n_events),   450, tol = 1)
chk("Gang Violence Repetition cont beta = 1.19", rep_gv$mean_b, 1.19, tol = 0.01)
chk("Gang Violence Repetition cont t = 26.3",    rep_gv$tstat,  26.3, tol = 0.2)
chk("Glances init_pct = 0.4%",
    100 * mean(gl_sub$n_init / gl_sub$n_events),  0.4, tol = 0.1)
chk("Gang Violence init_pct = 46.4%",
    100 * mean(gv_sub$n_init / gv_sub$n_events), 46.4, tol = 0.2)
chk("PMs init_pct = 42.4%",
    100 * mean(res_df$n_init[res_df$cat == "PMs"] /
               res_df$n_events[res_df$cat == "PMs"]), 42.4, tol = 0.2)

# ── Table 2: Model fit summary (tab:ll) ───────────────────────────────────────
section("Table 2: Model fit summary (tab:ll)")

paper_ll <- data.frame(
  domain   = c("Classroom","Email","Gang Violence","Glances","PMs",
                "Phone (DoD)","SMS","Twitter","War"),
  delta_aic = c(752, 2554, 549, 376, 7415, 220, 1442, 243, 489),
  acc_sep   = c(0.430, 0.262, 0.140, 0.090, 0.877, 0.424, 0.938, 0.337, 0.986),
  acc_pool  = c(0.390, 0.240, 0.075, 0.090, 0.628, 0.323, 0.879, 0.292, 0.932),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(paper_ll))) {
  d  <- paper_ll$domain[i]
  sub <- res_df[res_df$cat == d, ]
  chk(paste0(d, " delta_aic"),
      mean(sub$delta_aic, na.rm = TRUE),  paper_ll$delta_aic[i], tol = 5)
  chk(paste0(d, " acc_sep"),
      mean(sub$acc_sep_cont, na.rm = TRUE), paper_ll$acc_sep[i],  tol = 0.002)
  chk(paste0(d, " acc_pool"),
      mean(sub$acc_pool, na.rm = TRUE),   paper_ll$acc_pool[i],  tol = 0.002)
}

# AIC preferred in 210/210
chk("AIC: sep preferred in 210/210 networks",
    sum(res_df$delta_aic > 0, na.rm = TRUE), 210, tol = 0)

# Ecology-weighted mean accuracy
means_sep  <- sapply(domain_order, function(d) mean(res_df$acc_sep_cont[res_df$cat==d], na.rm=TRUE))
means_pool <- sapply(domain_order, function(d) mean(res_df$acc_pool[res_df$cat==d], na.rm=TRUE))
chk("Ecology-weighted mean acc_sep = 0.498",  mean(means_sep),  0.498, tol = 0.002)
chk("Ecology-weighted mean acc_pool = 0.428", mean(means_pool), 0.428, tol = 0.002)
chk("Sep > pool in all 9 domains",
    sum(means_sep > means_pool), 9, tol = 0)

# AIC visibility claim (simple means, not ecology-weighted)
low_vis  <- c("Email","Glances","PMs","SMS")
high_vis <- c("Classroom","Gang Violence","Twitter","Phone (DoD)","War")
low_mean  <- mean(unlist(lapply(low_vis,  function(d) res_df$delta_aic[res_df$cat==d])), na.rm=TRUE)
high_mean <- mean(unlist(lapply(high_vis, function(d) res_df$delta_aic[res_df$cat==d])), na.rm=TRUE)
chk("Low-visibility mean delta_AIC ~ 1637",  low_mean,  1637, tol = 10)
chk("High-visibility mean delta_AIC ~ 546",  high_mean,  546, tol = 10)

# Proportional gain range (Gang Violence and PMs)
gain_gv <- (means_sep["Gang Violence"] - means_pool["Gang Violence"]) / means_pool["Gang Violence"] * 100
gain_pm <- (means_sep["PMs"]           - means_pool["PMs"])           / means_pool["PMs"]           * 100
chk("Gang Violence proportional gain ~86% (paper says 40-88%)", gain_gv, 86, tol = 2)
chk("PMs proportional gain ~40% (paper says 40-88%)",           gain_pm, 40, tol = 2)

# ── Table 3: Triadic closure (tab:init) ───────────────────────────────────────
section("Table 3: Triadic closure (tab:init)")

# Initiation NOTPSnd
notp_i <- list(
  Email          = list(b = 1.67, t = 27.9),
  "Gang Violence"= list(b = 0.70, t = 38.1),
  Glances        = list(b = 0.19, t =  2.8),
  PMs            = list(b = 1.57, t =  9.4),
  "Phone (DoD)"  = list(b = 0.58, t = NA),
  SMS            = list(b = 2.83, t =  4.1),
  Twitter        = list(b = 1.64, t = NA),
  War            = list(b = 0.08, t = NA),
  Classroom      = list(b = 1.61, t = NA)
)
for (d in names(notp_i)) {
  v <- agg(res_df[res_df$cat == d, ], "init_NOTPSnd", "init_NOTPSnd_z")
  chk(paste0(d, " init NOTPSnd beta"),  v$mean_b,  notp_i[[d]]$b, tol = 0.01)
  if (!is.na(notp_i[[d]]$t))
    chk(paste0(d, " init NOTPSnd t"),   v$tstat,   notp_i[[d]]$t, tol = 0.3)
}

# Initiation NITPSnd (just multi-network domains)
nitp_i <- list(
  Email          = list(b =  0.43, t =  8.7),
  "Gang Violence"= list(b =  0.44, t = 21.0),
  PMs            = list(b = -0.19, t = -1.3),
  SMS            = list(b = -0.12, t = -0.2)
)
for (d in names(nitp_i)) {
  v <- agg(res_df[res_df$cat == d, ], "init_NITPSnd", "init_NITPSnd_z")
  chk(paste0(d, " init NITPSnd beta"),  v$mean_b,  nitp_i[[d]]$b, tol = 0.01)
  chk(paste0(d, " init NITPSnd t"),     v$tstat,   nitp_i[[d]]$t, tol = 0.3)
}

# Continuation NOTPSnd
notp_c <- list(
  Email          = list(b = 1.05, t = 19.9),
  "Gang Violence"= list(b = 0.70, t = 30.0),
  Glances        = list(b = 0.30, t =  5.1),
  PMs            = list(b = 0.30, t =  5.1),
  "Phone (DoD)"  = list(b = 0.97, t = NA),
  SMS            = list(b = 0.26, t =  1.6),
  Twitter        = list(b = 1.32, t = NA),
  War            = list(b = -0.22, t = NA),
  Classroom      = list(b = -0.04, t = NA)
)
for (d in names(notp_c)) {
  v <- agg(res_df[res_df$cat == d, ], "cont_NOTPSnd", "cont_NOTPSnd_z")
  chk(paste0(d, " cont NOTPSnd beta"),  v$mean_b,  notp_c[[d]]$b, tol = 0.01)
}

# NOTP init/cont ratios cited in text
ratios <- list(
  "Gang Violence" = list(r = 0.99, tol = 0.02),
  Email           = list(r = 1.59, tol = 0.03),
  PMs             = list(r = 5.25, tol = 0.10),
  SMS             = list(r = 10.94, tol = 0.20)
)
for (d in names(ratios)) {
  vi <- agg(res_df[res_df$cat == d, ], "init_NOTPSnd", "init_NOTPSnd_z")$mean_b
  vc <- agg(res_df[res_df$cat == d, ], "cont_NOTPSnd", "cont_NOTPSnd_z")$mean_b
  chk(paste0(d, " NOTP init/cont ratio"),
      vi / vc, ratios[[d]]$r, tol = ratios[[d]]$tol)
}

# ── Text: Reciprocity section ─────────────────────────────────────────────────
section("Text: Reciprocity — RecInitiation (initiation model)")

rec_init <- list(
  Email          = list(b =  1.815, t = 15.9),
  "Gang Violence"= list(b =  0.860, t = 11.2),
  PMs            = list(b =  2.748, t =  3.6),
  Glances        = list(b =  0.369, t =  3.1),
  SMS            = list(b = 14.157, t =  3.1),
  War            = list(b =  2.111, t =  1.9),
  Classroom      = list(b = -2.176, t = -1.8),
  Twitter        = list(b = -0.939, t = -0.4),
  "Phone (DoD)"  = list(b = -0.047, t = -0.03)
)
for (d in names(rec_init)) {
  v <- agg(res_df[res_df$cat == d, ], "init_RRecInitiation", "init_RRecInitiation_z")
  chk(paste0(d, " RecInitiation init beta"), v$mean_b, rec_init[[d]]$b, tol = 0.01)
  chk(paste0(d, " RecInitiation init t"),    v$tstat,  rec_init[[d]]$t, tol = 0.2)
}

section("Text: Reciprocity — sign flip rates (init vs pool RecInitiation)")

flip_paper <- list(
  Glances        = 0.60,
  Email          = 1.00,
  War            = 1.00,
  SMS            = 0.33,
  PMs            = 0.75,
  "Gang Violence"= 0.27
)
for (d in names(flip_paper)) {
  sub <- res_df[res_df$cat == d, ]
  ok  <- is.finite(sub$init_RRecInitiation) & is.finite(sub$pool_RRecInitiation)
  flip_rate <- mean(sign(sub$init_RRecInitiation[ok]) != sign(sub$pool_RRecInitiation[ok]))
  chk(paste0(d, " sign-flip rate"), flip_rate, flip_paper[[d]], tol = 0.02)
}

section("Text: RecRecencyLog continuation")

rrl_c <- list(
  Email          = list(b = -0.140, t = -12.1),
  "Gang Violence"= list(b = -0.236, t = -56.9),
  PMs            = list(b = -0.206, t =  -2.8),
  SMS            = list(b = -0.562, t =  -5.5),
  "Phone (DoD)"  = list(b = -0.292, t =  -5.4)
)
for (d in names(rrl_c)) {
  v <- agg(res_df[res_df$cat == d, ], "cont_RecRecencyLog", "cont_RecRecencyLog_z")
  chk(paste0(d, " RecRecencyLog cont beta"), v$mean_b, rrl_c[[d]]$b, tol = 0.005)
  chk(paste0(d, " RecRecencyLog cont t"),    v$tstat,  rrl_c[[d]]$t, tol = 0.5)
}

# ── Text: Continuation model ──────────────────────────────────────────────────
section("Text: Continuation model — Repetition")

rep_c <- list(
  Email          = list(b = 2.093, t = 16.8),
  "Gang Violence"= list(b = 1.192, t = 26.3),
  Glances        = list(b = 3.648, t = 51.8),
  PMs            = list(b = 0.343, t =  3.3),
  SMS            = list(b = 0.829, t =  8.6)
)
for (d in names(rep_c)) {
  v <- agg(res_df[res_df$cat == d, ], "cont_Repetition", "cont_Repetition_z")
  chk(paste0(d, " Repetition cont beta"), v$mean_b, rep_c[[d]]$b, tol = 0.01)
  chk(paste0(d, " Repetition cont t"),    v$tstat,  rep_c[[d]]$t, tol = 0.2)
}

section("Text: Continuation model — Imitation")
im_c <- list(
  Email   = list(b = -0.428, t = -9.7),
  Glances = list(b =  0.163, t =  5.0)
)
for (d in names(im_c)) {
  v <- agg(res_df[res_df$cat == d, ], "cont_Imitation", "cont_Imitation_z")
  chk(paste0(d, " Imitation cont beta"), v$mean_b, im_c[[d]]$b, tol = 0.005)
  chk(paste0(d, " Imitation cont t"),    v$tstat,  im_c[[d]]$t, tol = 0.3)
}

# ── Greedy forward selection ──────────────────────────────────────────────────
section("Text: Greedy forward selection — first-step statistics")

# Reconstruct cats and step-1 selections from greedy_sep (applies stopping rule)
cats <- sapply(greedy_sep, `[[`, "cat")
n_nets <- length(greedy_sep)

# Stopping-rule application: keep only steps where delta_aic > 0
# (greedy_sep stores full paths; stopping rule means we take the path
#  up to the last term that improved AIC — but the stored paths already
#  reflect the run order, not post-hoc stopping. The paper's stopping-rule
#  claim was verified in greedy_acc2.R. Here we use the step-1 selection
#  which is always the best single term regardless of stopping.)

step1_init <- sapply(greedy_sep, function(el) {
  p <- el$path_init
  if (is.null(p) || length(p) == 0) return(NA_character_)
  p[1]
})
step1_cont <- sapply(greedy_sep, function(el) {
  p <- el$path_cont
  if (is.null(p) || length(p) == 0) return(NA_character_)
  p[1]
})

# Overall first-step frequency
ri_frac_overall <- mean(grepl("RecInitiation", step1_init, fixed = TRUE), na.rm = TRUE)
sr_frac_overall <- mean(grepl("SndRecencyLog", step1_cont,  fixed = TRUE) |
                        grepl("SndRecencylog", step1_cont,  fixed = TRUE), na.rm = TRUE)
io_cont_frac    <- mean(grepl("IOviolation",   step1_cont,  fixed = TRUE), na.rm = TRUE)

chk("Init: RecInitiation first-selected overall = 46%", ri_frac_overall, 0.46, tol = 0.03)
chk("Cont: SndRecencylog first-selected overall = 56%", sr_frac_overall, 0.56, tol = 0.03)
chk("Cont: IOviolation first-selected overall = 33%",   io_cont_frac,    0.33, tol = 0.03)

# "7 of 9 domains" — RecInitiation leads init, SndRecencylog leads cont
init_dom_winner <- sapply(domain_order, function(d) {
  idx <- which(cats == d)
  if (length(idx) == 0) return(NA_character_)
  s1 <- step1_init[idx]
  s1 <- s1[!is.na(s1)]
  if (length(s1) == 0) return(NA_character_)
  names(sort(table(s1), decreasing = TRUE))[1]
})
cont_dom_winner <- sapply(domain_order, function(d) {
  idx <- which(cats == d)
  if (length(idx) == 0) return(NA_character_)
  s1 <- step1_cont[idx]
  s1 <- s1[!is.na(s1)]
  if (length(s1) == 0) return(NA_character_)
  names(sort(table(s1), decreasing = TRUE))[1]
})

n_ri_dom <- sum(grepl("RecInitiation", init_dom_winner, fixed = TRUE), na.rm = TRUE)
n_sr_dom <- sum(grepl("SndRecencyLog|SndRecencylog", cont_dom_winner), na.rm = TRUE)
chk("Init: RecInitiation leads in 7 of 9 domains", n_ri_dom, 7, tol = 0)
chk("Cont: SndRecencylog leads in 7 of 9 domains", n_sr_dom, 7, tol = 0)

# Domain-level first-step % for RecInitiation (cited in text)
dom_ri_pcts <- list(
  Email          = 0.52,
  "Gang Violence"= 0.90,
  PMs            = 0.75,
  "Phone (DoD)"  = 1.00,
  SMS            = 0.67,
  War            = 1.00
)
for (d in names(dom_ri_pcts)) {
  idx <- which(cats == d)
  if (length(idx) == 0) { cat(sprintf("  %-65s  [NO DATA]\n", paste0(d, " RecInit step-1 %"))); next }
  s1 <- step1_init[idx]
  frac <- mean(grepl("RecInitiation", s1, fixed = TRUE), na.rm = TRUE)
  chk(paste0(d, " RecInit step-1 % = ", round(100 * dom_ri_pcts[[d]]), "%"),
      frac, dom_ri_pcts[[d]], tol = 0.05)
}

# NOTPSnd in first 3 steps: 43% of init fits
notp_in_3_init <- mean(sapply(greedy_sep, function(el) {
  p <- el$path_init; if (is.null(p)) return(NA)
  any(grepl("NOTPSnd", head(p, 3), fixed = TRUE))
}), na.rm = TRUE)
notp_in_3_cont <- mean(sapply(greedy_sep, function(el) {
  p <- el$path_cont; if (is.null(p)) return(NA)
  any(grepl("NOTPSnd", head(p, 3), fixed = TRUE))
}), na.rm = TRUE)
chk("Init: NOTPSnd in first 3 steps in 43% of fits",  notp_in_3_init, 0.43, tol = 0.04)
chk("Cont: NOTPSnd in first 3 steps in <3% of fits",  notp_in_3_cont, 0.03, tol = 0.03)

# ── Appendix B: Ecology-weighted table (tab:weighted) ────────────────────────
section("Appendix B: Ecology-weighted table (tab:weighted)")

paper_wt <- data.frame(
  domain    = c("Classroom","Email","Gang Violence","Glances","PMs",
                "Phone (DoD)","SMS","Twitter","War"),
  pref      = rep(1.00, 9),
  aic_per_ev = c(1.088, 0.498, 1.225, 0.038, 1.064, 1.106, 0.386, 1.039, 0.279),
  acc_sep   = c(0.430, 0.263, 0.139, 0.090, 0.875, 0.389, 0.936, 0.376, 0.979),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(paper_wt))) {
  d   <- paper_wt$domain[i]
  sub <- res_df[res_df$cat == d, ]
  pref_rate  <- mean(sub$delta_aic > 0, na.rm = TRUE)
  aic_per_ev <- mean(sub$delta_aic / sub$n_events, na.rm = TRUE)
  acc_sep    <- mean(sub$acc_sep_cont, na.rm = TRUE)
  chk(paste0(d, " pref = 1.00"),        pref_rate,  1.00,                   tol = 0)
  chk(paste0(d, " AIC/ev"),             aic_per_ev, paper_wt$aic_per_ev[i], tol = 0.005)
  chk(paste0(d, " acc_sep (tab:wt)"),   acc_sep,    paper_wt$acc_sep[i],    tol = 0.002)
}

wt_overall_sep <- mean(sapply(domain_order, function(d)
  mean(res_df$acc_sep_cont[res_df$cat == d], na.rm = TRUE)))
chk("tab:weighted Wtd overall acc_sep = 0.498", wt_overall_sep, 0.498, tol = 0.002)

# ── Appendix C: Independent vs Interacted (tab:acc) ──────────────────────────
section("Appendix C: tab:acc — accuracy by specification")

paper_acc <- data.frame(
  domain   = c("Classroom","Email","Gang Violence","Glances","PMs",
               "Phone (DoD)","SMS","Twitter","War"),
  acc_sep  = c(0.430, 0.262, 0.140, 0.090, 0.877, 0.424, 0.938, 0.337, 0.986),
  acc_pool = c(0.390, 0.240, 0.075, 0.090, 0.628, 0.323, 0.879, 0.292, 0.932),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(paper_acc))) {
  d   <- paper_acc$domain[i]
  sub <- res_df[res_df$cat == d, ]
  chk(paste0(d, " acc_sep (tab:acc)"),
      mean(sub$acc_sep_cont, na.rm = TRUE), paper_acc$acc_sep[i],  tol = 0.002)
  chk(paste0(d, " acc_pool (tab:acc)"),
      mean(sub$acc_pool,     na.rm = TRUE), paper_acc$acc_pool[i], tol = 0.002)
}

# Interacted domains: converged means (Email, War, SMS, PMs)
conv_dom <- c("Email","War","SMS","PMs")
sep_conv  <- mean(sapply(conv_dom, function(d)
  mean(res_df$acc_sep_cont[res_df$cat == d], na.rm = TRUE)))
chk("Converged domains ecology-wtd sep = 0.766", sep_conv, 0.766, tol = 0.003)

# ── Appendix D: Coefficient tables ───────────────────────────────────────────
section("Appendix D: Coefficient tables (coef_table_sep.tex spot-checks)")

# These match the auto-generated table values from coef_table.R
coef_checks <- list(
  list(d="Glances",       m="init", bc="init_RRecInitiation", zc="init_RRecInitiation_z", b=0.369, stars=2L),
  list(d="Email",         m="init", bc="init_RRecInitiation", zc="init_RRecInitiation_z", b=1.815, stars=3L),
  list(d="SMS",           m="init", bc="init_RRecInitiation", zc="init_RRecInitiation_z", b=14.157, stars=2L),
  list(d="PMs",           m="init", bc="init_RRecInitiation", zc="init_RRecInitiation_z", b=2.748, stars=3L),
  list(d="Gang Violence", m="init", bc="init_RRecInitiation", zc="init_RRecInitiation_z", b=0.860, stars=3L),
  list(d="Glances",       m="cont", bc="cont_Repetition",     zc="cont_Repetition_z",     b=3.648, stars=3L),
  list(d="Email",         m="cont", bc="cont_Repetition",     zc="cont_Repetition_z",     b=2.093, stars=3L),
  list(d="Gang Violence", m="cont", bc="cont_Repetition",     zc="cont_Repetition_z",     b=1.192, stars=3L),
  list(d="SMS",           m="cont", bc="cont_IOviolation",    zc="cont_IOviolation_z",    b=-4.143, stars=3L),
  list(d="PMs",           m="cont", bc="cont_IOviolation",    zc="cont_IOviolation_z",    b=-3.494, stars=3L)
)
for (cc in coef_checks) {
  v <- agg(res_df[res_df$cat == cc$d, ], cc$bc, cc$zc)
  chk(paste0(cc$d, " ", cc$m, " ", sub(".*_([^_]+)$","\\1",cc$bc), " beta"),
      v$mean_b, cc$b, tol = 0.005)
}

# ── Greedy accuracy equivalence (Limitations paragraph) ──────────────────────
section("Limitations: greedy vs fixed-spec accuracy equivalence")

# greedy_acc2_results.Rdata: compare fixed-spec and greedy accuracy on 78-network subsample
# The paper claims identical to 4 decimal places
ga2_domains <- sort(unique(ga2$cat))
for (d in ga2_domains) {
  sub <- ga2[ga2$cat == d, ]
  fixed_mean  <- mean(sub$acc_sep_cont,    na.rm = TRUE)
  greedy_mean <- mean(sub$acc_greedy_cont, na.rm = TRUE)
  diff <- abs(fixed_mean - greedy_mean)
  ok <- diff < 5e-4
  cat(sprintf("  %-20s  fixed=%.4f  greedy=%.4f  diff=%.5f  %s\n",
    d, fixed_mean, greedy_mean, diff,
    if (ok) "OK (<0.0005)" else "MISMATCH"))
}

# ── Summary ───────────────────────────────────────────────────────────────────
section("Done")
cat("Scroll up for any MISMATCH lines.\n")
cat("tol for beta comparisons: 0.005–0.01; for t-stats: 0.2–0.5.\n")
cat("Larger deviations indicate the paper text or table needs correction.\n")
