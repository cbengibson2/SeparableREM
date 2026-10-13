# coef_table.R
#
# Produces two LaTeX tables of ecology-weighted mean coefficients with
# significance stars:
#   1. Separable model: initiation and continuation sub-models side by side
#   2. Pooled model
#
# Reads from revision_results.Rdata (the authoritative run).
#
# Multi-network domains: t-stat = mean_beta / (sd_beta / sqrt(N))
# Single-network domains: z-stat from the fitted model summary
#
# Output:
#   paper/coef_table_sep.tex
#   paper/coef_table_pool.tex

rm(list = ls())

load("revision_results.Rdata")

domain_order <- c("Glances", "Email", "War", "SMS", "Classroom",
                  "Twitter", "Phone (DoD)", "PMs", "Gang Violence")

# ── Aggregation ───────────────────────────────────────────────────────────────

agg_domain <- function(df_sub, beta_col, z_col) {
  bvals <- df_sub[[beta_col]]
  zvals <- df_sub[[z_col]]
  bvals <- bvals[is.finite(bvals)]
  zvals <- zvals[is.finite(zvals)]
  N <- length(bvals)
  if (N == 0) return(list(mean_b = NA_real_, tstat = NA_real_, N = 0L))
  mean_b <- mean(bvals)
  if (N == 1) {
    tstat <- if (length(zvals) > 0) zvals[1] else NA_real_
  } else {
    sd_b  <- sd(bvals)
    tstat <- if (!is.na(sd_b) && sd_b > 0) mean_b / (sd_b / sqrt(N)) else NA_real_
  }
  list(mean_b = mean_b, tstat = tstat, N = N)
}

# ── Cell renderer ─────────────────────────────────────────────────────────────

stars <- function(t) {
  if (is.na(t)) return("")
  at <- abs(t)
  if (at > 3.29) return("^{***}")
  if (at > 2.58) return("^{**}")
  if (at > 1.96) return("^{*}")
  return("")
}

cell <- function(mean_b, tstat) {
  if (is.na(mean_b)) return("\\multicolumn{1}{c}{---}")
  st <- stars(tstat)
  sprintf("$%s%.3f%s$",
          if (mean_b >= 0) "\\phantom{-}" else "",
          mean_b, st)
}

# ── Term definitions ──────────────────────────────────────────────────────────

INIT_TERMS <- list(
  RRecInit  = c("init_RRecInitiation", "init_RRecInitiation_z"),
  NOTP      = c("init_NOTPSnd",        "init_NOTPSnd_z"),
  NITP      = c("init_NITPSnd",        "init_NITPSnd_z"),
  RecRec    = c("init_RecRecencyLog",  "init_RecRecencyLog_z"),
  Imit      = c("init_Imitation",      "init_Imitation_z"),
  IOviol    = c("init_IOviolation",    "init_IOviolation_z")
)

CONT_TERMS <- list(
  Rep       = c("cont_Repetition",    "cont_Repetition_z"),
  NOTP      = c("cont_NOTPSnd",       "cont_NOTPSnd_z"),
  NITP      = c("cont_NITPSnd",       "cont_NITPSnd_z"),
  RecRec    = c("cont_RecRecencyLog", "cont_RecRecencyLog_z"),
  Imit      = c("cont_Imitation",     "cont_Imitation_z"),
  IOviol    = c("cont_IOviolation",   "cont_IOviolation_z")
)

POOL_TERMS <- list(
  SndInit   = c("pool_SndInitiation",  "pool_SndInitiation_z"),
  RRecInit  = c("pool_RRecInitiation", "pool_RRecInitiation_z"),
  Rep       = c("pool_Repetition",     "pool_Repetition_z"),
  NOTP      = c("pool_NOTPSnd",        "pool_NOTPSnd_z"),
  NITP      = c("pool_NITPSnd",        "pool_NITPSnd_z"),
  RecRec    = c("pool_RecRecencyLog",  "pool_RecRecencyLog_z"),
  Imit      = c("pool_Imitation",      "pool_Imitation_z"),
  IOviol    = c("pool_IOviolation",    "pool_IOviolation_z")
)

n_i <- length(INIT_TERMS)
n_c <- length(CONT_TERMS)
n_p <- length(POOL_TERMS)

cap_note_shared <- paste0(
  "Ecology-weighted mean coefficients. ",
  "For multi-network domains, significance stars reflect the ecology-weighted ",
  "$t$-statistic ($\\bar{\\beta}/(\\hat{\\sigma}/\\sqrt{N})$); ",
  "for single-network domains ($^{\\ddagger}$) the individual model $z$-statistic is used. ",
  "$^{*}p<.05$, $^{**}p<.01$, $^{***}p<.001$. ",
  "\\textit{RRecInit}: RecInitiation; \\textit{Rep}: Repetition; ",
  "\\textit{NOTP}: NOTPSnd; \\textit{NITP}: NITPSnd; ",
  "\\textit{RecRec}: RecRecencyLog; \\textit{Imit}: Imitation; ",
  "\\textit{IOviol}: IOviolation; \\textit{SndInit}: SndInitiation. ",
  "Gang Violence ($^{\\dagger}$) is synthetically generated."
)

# ── Table 1: Separable — split into initiation panel and continuation panel ───

sup_label <- function(dom) {
  label <- if (dom == "Private Messages") "PMs" else dom
  sup   <- if (dom %in% c("Classroom","Phone (DoD)","Twitter","War")) "$^{\\ddagger}$"
           else if (dom == "Gang Violence") "$^{\\dagger}$"
           else ""
  list(label = label, sup = sup)
}

make_row_init <- function(dom) {
  sub  <- res_df[res_df$cat == dom, ]
  N    <- nrow(sub)
  ls   <- sup_label(dom)
  cells <- sapply(INIT_TERMS, function(cols)
    cell(agg_domain(sub, cols[1], cols[2])$mean_b,
         agg_domain(sub, cols[1], cols[2])$tstat))
  paste0(sprintf("%s%s & %3d", ls$label, ls$sup, N), " & ",
         paste(cells, collapse = " & "), " \\\\")
}

make_row_cont <- function(dom) {
  sub  <- res_df[res_df$cat == dom, ]
  N    <- nrow(sub)
  ls   <- sup_label(dom)
  cells <- sapply(CONT_TERMS, function(cols)
    cell(agg_domain(sub, cols[1], cols[2])$mean_b,
         agg_domain(sub, cols[1], cols[2])$tstat))
  paste0(sprintf("%s%s & %3d", ls$label, ls$sup, N), " & ",
         paste(cells, collapse = " & "), " \\\\")
}

col_spec_i <- paste0("lr", paste(rep("r", n_i), collapse = ""))
col_spec_c <- paste0("lr", paste(rep("r", n_c), collapse = ""))
hdr_i <- paste(names(INIT_TERMS), collapse = " & ")
hdr_c <- paste(names(CONT_TERMS), collapse = " & ")

tex_init <- c(
  "% Generated by coef_table.R (separable model) — do not edit by hand",
  "\\begin{table}[p]",
  "\\centering",
  sprintf("\\caption{Separable REM coefficient estimates by domain: initiation sub-model. %s}",
          cap_note_shared),
  "\\label{tab:coef_sep}",
  "{\\small",
  paste0("\\begin{tabular}{", col_spec_i, "}"),
  "\\toprule",
  paste0("\\textbf{Domain} & $N$ & ", hdr_i, " \\\\"),
  "\\midrule",
  paste(sapply(domain_order, make_row_init), collapse = "\n"),
  "\\bottomrule",
  "\\end{tabular}",
  "}",
  "\\end{table}"
)

tex_cont <- c(
  "\\begin{table}[p]",
  "\\centering",
  sprintf("\\caption{Separable REM coefficient estimates by domain: continuation sub-model. %s}",
          cap_note_shared),
  "\\label{tab:coef_cont}",
  "{\\small",
  paste0("\\begin{tabular}{", col_spec_c, "}"),
  "\\toprule",
  paste0("\\textbf{Domain} & $N$ & ", hdr_c, " \\\\"),
  "\\midrule",
  paste(sapply(domain_order, make_row_cont), collapse = "\n"),
  "\\bottomrule",
  "\\end{tabular}",
  "}",
  "\\end{table}"
)

writeLines(c(tex_init, "", tex_cont),
           "paper/coef_table_sep.tex")
cat("Written: paper/coef_table_sep.tex\n")

# ── Table 2: Pooled model ─────────────────────────────────────────────────────

make_row_pool <- function(dom) {
  sub   <- res_df[res_df$cat == dom, ]
  N     <- nrow(sub)
  label <- if (dom == "Private Messages") "PMs" else dom
  sup   <- if (dom %in% c("Classroom","Phone (DoD)","Twitter","War")) "$^{\\ddagger}$"
           else if (dom == "Gang Violence") "$^{\\dagger}$"
           else ""

  cells <- sapply(POOL_TERMS, function(cols)
    cell(agg_domain(sub, cols[1], cols[2])$mean_b,
         agg_domain(sub, cols[1], cols[2])$tstat))

  paste0(sprintf("%s%s & %3d", label, sup, N), " & ",
         paste(cells, collapse = " & "), " \\\\")
}

col_spec_pool <- paste0("lr", paste(rep("r", n_p), collapse = ""))
hdr_p <- paste(names(POOL_TERMS), collapse = " & ")

tex_pool <- c(
  "% Generated by coef_table.R (pooled model) — do not edit by hand",
  "\\begin{table}[p]",
  "\\centering",
  sprintf("\\caption{Pooled REM coefficient estimates by domain. %s}", cap_note_shared),
  "\\label{tab:coef_pool}",
  "{\\small",
  paste0("\\begin{tabular}{", col_spec_pool, "}"),
  "\\toprule",
  paste0("\\textbf{Domain} & $N$ & ", hdr_p, " \\\\"),
  "\\midrule",
  paste(sapply(domain_order, make_row_pool), collapse = "\n"),
  "\\bottomrule",
  "\\end{tabular}",
  "}",
  "\\end{table}"
)

writeLines(tex_pool, "paper/coef_table_pool.tex")
cat("Written: paper/coef_table_pool.tex\n")

# ── Console preview ───────────────────────────────────────────────────────────

cat("\n=== Separable model — Initiation ===\n")
cat(sprintf("%-20s %3s  %s\n", "Domain", "N",
    paste(names(INIT_TERMS), collapse="       ")))
for (dom in domain_order) {
  sub <- res_df[res_df$cat==dom,]; N <- nrow(sub)
  vals <- sapply(INIT_TERMS, function(cols) {
    ag <- agg_domain(sub, cols[1], cols[2])
    if (is.na(ag$mean_b)) "   ---  " else
      sprintf("%6.3f%s", ag$mean_b, if (!is.na(ag$tstat) && abs(ag$tstat)>1.96) "*" else " ")
  })
  cat(sprintf("%-20s %3d  %s\n", dom, N, paste(vals, collapse="  ")))
}

cat("\n=== Separable model — Continuation ===\n")
cat(sprintf("%-20s %3s  %s\n", "Domain", "N",
    paste(names(CONT_TERMS), collapse="       ")))
for (dom in domain_order) {
  sub <- res_df[res_df$cat==dom,]; N <- nrow(sub)
  vals <- sapply(CONT_TERMS, function(cols) {
    ag <- agg_domain(sub, cols[1], cols[2])
    if (is.na(ag$mean_b)) "   ---  " else
      sprintf("%6.3f%s", ag$mean_b, if (!is.na(ag$tstat) && abs(ag$tstat)>1.96) "*" else " ")
  })
  cat(sprintf("%-20s %3d  %s\n", dom, N, paste(vals, collapse="  ")))
}

cat("\n=== Pooled model ===\n")
cat(sprintf("%-20s %3s  %s\n", "Domain", "N",
    paste(names(POOL_TERMS), collapse="       ")))
for (dom in domain_order) {
  sub <- res_df[res_df$cat==dom,]; N <- nrow(sub)
  vals <- sapply(POOL_TERMS, function(cols) {
    ag <- agg_domain(sub, cols[1], cols[2])
    if (is.na(ag$mean_b)) "   ---  " else
      sprintf("%6.3f%s", ag$mean_b, if (!is.na(ag$tstat) && abs(ag$tstat)>1.96) "*" else " ")
  })
  cat(sprintf("%-20s %3d  %s\n", dom, N, paste(vals, collapse="  ")))
}

cat("\nDone.\n")
