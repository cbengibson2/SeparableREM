# sign_table.R
#
# Produces a three-panel sign+significance table (LaTeX) for all statistics
# across all three models (initiation, continuation, pooled) by domain.
#
# Reads coefficients and z-stats from revision_results.Rdata (the authoritative
# run) rather than re-fitting, so the table is guaranteed to match the paper text.
#
# For multi-network domains: ecology-weighted t-stat = mean_beta / (sd/sqrt(N))
# For single-network domains: z-stat from the fitted model summary
#
# Output: paper/sign_table.tex (standalone tabular, ready to \input)

rm(list = ls())

load("revision_results.Rdata")

domain_order <- c("Glances", "Email", "War", "SMS", "Classroom",
                  "Twitter", "Phone (DoD)", "PMs", "Gang Violence")

# Remap "Private Messages" -> "PMs" for display (res_df uses "PMs" already)
cat(sprintf("Loaded %d networks across %d domains\n", nrow(res_df),
            length(unique(res_df$cat))))

# ── Column mapping: display name -> (beta_col, z_col) ─────────────────────────

INIT_TERMS <- list(
  RRecInitiation = c("init_RRecInitiation",   "init_RRecInitiation_z"),
  NOTPSnd        = c("init_NOTPSnd",           "init_NOTPSnd_z"),
  NITPSnd        = c("init_NITPSnd",           "init_NITPSnd_z"),
  RecRecencyLog  = c("init_RecRecencyLog",     "init_RecRecencyLog_z"),
  Imitation      = c("init_Imitation",         "init_Imitation_z"),
  IOviolation    = c("init_IOviolation",        "init_IOviolation_z")
)

CONT_TERMS <- list(
  Repetition    = c("cont_Repetition",    "cont_Repetition_z"),
  NOTPSnd       = c("cont_NOTPSnd",       "cont_NOTPSnd_z"),
  NITPSnd       = c("cont_NITPSnd",       "cont_NITPSnd_z"),
  RecRecencyLog = c("cont_RecRecencyLog", "cont_RecRecencyLog_z"),
  Imitation     = c("cont_Imitation",     "cont_Imitation_z"),
  IOviolation   = c("cont_IOviolation",   "cont_IOviolation_z")
)

POOL_TERMS <- list(
  SndInitiation  = c("pool_SndInitiation",  "pool_SndInitiation_z"),
  RRecInitiation = c("pool_RRecInitiation", "pool_RRecInitiation_z"),
  Repetition     = c("pool_Repetition",     "pool_Repetition_z"),
  NOTPSnd        = c("pool_NOTPSnd",        "pool_NOTPSnd_z"),
  NITPSnd        = c("pool_NITPSnd",        "pool_NITPSnd_z"),
  RecRecencyLog  = c("pool_RecRecencyLog",  "pool_RecRecencyLog_z"),
  Imitation      = c("pool_Imitation",      "pool_Imitation_z"),
  IOviolation    = c("pool_IOviolation",    "pool_IOviolation_z")
)

# ── Aggregation ───────────────────────────────────────────────────────────────
#
# Multi-network domain (N>1): t-stat = mean_beta / (sd_beta / sqrt(N))
# Single-network domain (N=1): use the model z-stat directly

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
  if (at > 3.29) return("***")
  if (at > 2.58) return("**")
  if (at > 1.96) return("*")
  return("")
}

cell <- function(mean_b, tstat) {
  if (is.na(mean_b)) return("---")
  if (abs(mean_b) < 0.02) return("$\\cdot$")
  sgn <- if (mean_b > 0) "$+$" else "$-$"
  st  <- stars(tstat)
  if (nchar(st) > 0) paste0(sgn, "\\textsuperscript{", st, "}") else sgn
}

# ── Build table rows ──────────────────────────────────────────────────────────

make_row <- function(cat_label, init_terms, cont_terms, pool_terms) {
  sub   <- res_df[res_df$cat == cat_label, ]
  N     <- nrow(sub)
  label <- if (cat_label == "Private Messages") "PMs" else cat_label

  cells_i <- sapply(names(init_terms), function(tm) {
    cols <- init_terms[[tm]]
    ag   <- agg_domain(sub, cols[1], cols[2])
    cell(ag$mean_b, ag$tstat)
  })
  cells_c <- sapply(names(cont_terms), function(tm) {
    cols <- cont_terms[[tm]]
    ag   <- agg_domain(sub, cols[1], cols[2])
    cell(ag$mean_b, ag$tstat)
  })
  cells_p <- sapply(names(pool_terms), function(tm) {
    cols <- pool_terms[[tm]]
    ag   <- agg_domain(sub, cols[1], cols[2])
    cell(ag$mean_b, ag$tstat)
  })

  paste0(
    sprintf("%-20s & %3d", label, N), " & ",
    paste(cells_i, collapse = " & "), " & ",
    paste(cells_c, collapse = " & "), " & ",
    paste(cells_p, collapse = " & "), " \\\\"
  )
}

# ── Column header abbreviations ───────────────────────────────────────────────

abbrev_i <- c("RRecInit", "NOTP", "NITP", "RecRec", "Imit", "IOviol")
abbrev_c <- c("Rep",      "NOTP", "NITP", "RecRec", "Imit", "IOviol")
abbrev_p <- c("SndInit", "RRecInit", "Rep", "NOTP", "NITP", "RecRec", "Imit", "IOviol")

n_i <- length(INIT_TERMS)
n_c <- length(CONT_TERMS)
n_p <- length(POOL_TERMS)

# ── Assemble LaTeX ────────────────────────────────────────────────────────────
#
# Split into two tables:
#   Panel A: Initiation + Continuation (14 columns)
#   Panel B: Pooled (10 columns)

cap_note <- paste0(
  "Sign and significance of ecology-weighted coefficient estimates by domain. ",
  "For multi-network domains the reported sign is that of the domain mean; significance stars reflect ",
  "the ecology-weighted $t$-statistic ($\\bar{\\beta}/(\\hat{\\sigma}/\\sqrt{N})$). ",
  "For single-network domains the sign and stars are taken from the individual model $z$-statistic. ",
  "$+$/$-$: positive/negative; $\\cdot$: near zero ($|\\bar{\\beta}| < 0.02$); ---: term not identified. ",
  "Significance: ${}^{*}p<.05$, ${}^{**}p<.01$, ${}^{***}p<.001$. ",
  "\\textit{RRecInit}: RecInitiation; \\textit{Rep}: Repetition; \\textit{NOTP}: NOTPSnd; ",
  "\\textit{NITP}: NITPSnd; \\textit{RecRec}: RecRecencylog; \\textit{Imit}: Imitation; ",
  "\\textit{IOviol}: IOviolation; \\textit{SndInit}: SndInitiation."
)

make_row_ic <- function(cat_label) {
  sub   <- res_df[res_df$cat == cat_label, ]
  N     <- nrow(sub)
  label <- if (cat_label == "Private Messages") "PMs" else cat_label
  cells_i <- sapply(names(INIT_TERMS), function(tm) {
    cols <- INIT_TERMS[[tm]]; ag <- agg_domain(sub, cols[1], cols[2])
    cell(ag$mean_b, ag$tstat)
  })
  cells_c <- sapply(names(CONT_TERMS), function(tm) {
    cols <- CONT_TERMS[[tm]]; ag <- agg_domain(sub, cols[1], cols[2])
    cell(ag$mean_b, ag$tstat)
  })
  paste0(sprintf("%-20s & %3d", label, N), " & ",
         paste(cells_i, collapse = " & "), " & ",
         paste(cells_c, collapse = " & "), " \\\\")
}

make_row_p <- function(cat_label) {
  sub   <- res_df[res_df$cat == cat_label, ]
  N     <- nrow(sub)
  label <- if (cat_label == "Private Messages") "PMs" else cat_label
  cells_p <- sapply(names(POOL_TERMS), function(tm) {
    cols <- POOL_TERMS[[tm]]; ag <- agg_domain(sub, cols[1], cols[2])
    cell(ag$mean_b, ag$tstat)
  })
  paste0(sprintf("%-20s & %3d", label, N), " & ",
         paste(cells_p, collapse = " & "), " \\\\")
}

header_i <- paste(abbrev_i, collapse = " & ")
header_c <- paste(abbrev_c, collapse = " & ")
header_p <- paste(abbrev_p, collapse = " & ")

# Panel A: Initiation + Continuation
col_spec_ic <- paste0("lr",
                      paste(rep("c", n_i), collapse = ""),
                      paste(rep("c", n_c), collapse = ""))
span_i  <- sprintf("\\multicolumn{%d}{c}{\\textbf{Initiation}}", n_i)
span_c  <- sprintf("\\multicolumn{%d}{c}{\\textbf{Continuation}}", n_c)
cmidr_i <- sprintf("\\cmidrule(lr){3-%d}", 2 + n_i)
cmidr_c <- sprintf("\\cmidrule(lr){%d-%d}", 3 + n_i, 2 + n_i + n_c)

tex_a <- c(
  "% Generated by sign_table.R — do not edit by hand",
  "\\begin{table}[p]",
  "\\centering",
  sprintf("\\caption{%s (Separable sub-models.)}", cap_note),
  "\\label{tab:signs}",
  "{\\small",
  paste0("\\begin{tabular}{", col_spec_ic, "}"),
  "\\toprule",
  paste0("\\textbf{Domain} & $N$ & ", span_i, " & ", span_c, " \\\\"),
  paste(cmidr_i, cmidr_c),
  paste0(" & & ", header_i, " & ", header_c, " \\\\"),
  "\\midrule",
  paste(sapply(domain_order, make_row_ic), collapse = "\n"),
  "\\bottomrule",
  "\\end{tabular}",
  "}",
  "\\end{table}"
)

# Panel B: Pooled
col_spec_p <- paste0("lr", paste(rep("c", n_p), collapse = ""))
span_p  <- sprintf("\\multicolumn{%d}{c}{\\textbf{Pooled}}", n_p)

tex_b <- c(
  "\\begin{table}[p]",
  "\\centering",
  sprintf("\\caption{%s (Pooled model, continued.)}", cap_note),
  "\\label{tab:signs_pool}",
  "{\\small",
  paste0("\\begin{tabular}{", col_spec_p, "}"),
  "\\toprule",
  paste0("\\textbf{Domain} & $N$ & ", header_p, " \\\\"),
  "\\midrule",
  paste(sapply(domain_order, make_row_p), collapse = "\n"),
  "\\bottomrule",
  "\\end{tabular}",
  "}",
  "\\end{table}"
)

out_path <- "paper/sign_table.tex"
writeLines(c(tex_a, "", tex_b), out_path)
cat(sprintf("\nLaTeX written to %s\n", out_path))

# ── Console summary ───────────────────────────────────────────────────────────

cat("\n=== Sign and significance by domain ===\n")
cat(sprintf("%-20s  %3s  %-20s  %-20s  %-20s\n",
            "Domain", "N", "init_RRecInit", "cont_Rep", "init_NOTP"))
for (dom in domain_order) {
  sub <- res_df[res_df$cat == dom, ]
  N   <- nrow(sub)
  ag_ri  <- agg_domain(sub, "init_RRecInitiation", "init_RRecInitiation_z")
  ag_rep <- agg_domain(sub, "cont_Repetition",     "cont_Repetition_z")
  ag_nt  <- agg_domain(sub, "init_NOTPSnd",        "init_NOTPSnd_z")
  cat(sprintf("%-20s  %3d  b=%6.3f t=%6.1f    b=%6.3f t=%6.1f    b=%6.3f t=%6.1f\n",
    dom, N,
    ifelse(is.na(ag_ri$mean_b),  0, ag_ri$mean_b),  ifelse(is.na(ag_ri$tstat),  0, ag_ri$tstat),
    ifelse(is.na(ag_rep$mean_b), 0, ag_rep$mean_b), ifelse(is.na(ag_rep$tstat), 0, ag_rep$tstat),
    ifelse(is.na(ag_nt$mean_b),  0, ag_nt$mean_b),  ifelse(is.na(ag_nt$tstat),  0, ag_nt$tstat)
  ))
}

cat("\nDone.\n")
