# greedy_opt_sep.R
#
# Separable greedy forward selection for REMs across all available datasets.
#
# For each network the event history is partitioned by SndInitiation:
#   Initiation subset  (SndInitiation == 0): greedy_forward on pe_init
#   Continuation subset (SndInitiation == 1): greedy_forward on pe_cont
#
# Each regime has its own candidate term pool. Terms that are identically
# zero in a regime (e.g. Repetition when SndInitiation == 0) are excluded
# from that regime's pool. The selection procedure and AIC criterion are
# identical to greedy_opt.R.

rm(list = ls())
devtools::load_all("llrem", quiet = TRUE)
library(parallel)
library(survival)
library(ggplot2)
library(tidyr)
library(dplyr)

# ── Paths ─────────────────────────────────────────────────────────────────────

DATA_DIR <- Sys.getenv("DATA_DIR", "")   # path to the DataP network data directory
OUT_DIR  <- Sys.getenv("OUT_DIR",  ".")  # output directory (default: project root)

# ── Domain metadata ───────────────────────────────────────────────────────────

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

# Subsample large domains for speed
set.seed(8214L)
gang_idx <- which(all_files$cat == "Gang Violence")
ftf_idx  <- which(all_files$cat == "Glances")
keep_idx <- c(
  setdiff(seq_len(nrow(all_files)), c(gang_idx, ftf_idx)),
  sample(gang_idx, min(10L, length(gang_idx))),
  sample(ftf_idx,  min(20L, length(ftf_idx)))
)
all_files <- all_files[sort(keep_idx), ]
cat(sprintf("After subsampling: %d files\n", nrow(all_files)))
print(table(all_files$cat))
cat("\n")

# ── Candidate term pools ──────────────────────────────────────────────────────
#
# Initiation pool: terms that are meaningful when i has never sent to j.
#   Repetition is excluded (identically zero).
#   SenderRecencyLog is excluded (log of last i→j time × SndInitiation = 0).
#   RRecInitiation is included (j may have sent to i even before first i→j).
#
# Continuation pool: terms meaningful when i has previously sent to j.
#   RRecInitiation excluded: near-constant for established dyads.
#   Repetition and SenderRecencyLog are continuation-specific.

term_groups_init <- list(
  RRecInitiation   = "RRecInitiation",
  Recencylog       = "Recencylog",      # j→i reverse recency; valid pre-first-contact
  WeightedRec      = "WeightedRec",
  Imitation        = "Imitation",
  IOviolation      = "IOviolation",
  NOSPSnd          = "NOSPSnd",
  NISPSnd          = "NISPSnd",
  NOTPSnd          = "NOTPSnd",
  NITPSnd          = "NITPSnd",
  OutDegSnd        = "OutDegSnd",
  InDegRec         = "InDegRec"
)

term_groups_cont <- list(
  Repetition       = "Repetition",
  WeightedRep      = "WeightedRep",
  DyadRecencyLog   = "DyadRecencyLog",  # -log1p(last i→j time)
  WeightedRec      = "WeightedRec",
  Recencylog       = "Recencylog",      # j→i reverse recency
  Imitation        = "Imitation",
  IOviolation      = "IOviolation",
  NOSPSnd          = "NOSPSnd",
  NISPSnd          = "NISPSnd",
  NOTPSnd          = "NOTPSnd",
  NITPSnd          = "NITPSnd",
  OutDegSnd        = "OutDegSnd",
  InDegRec         = "InDegRec"
)

# ── Greedy forward selection ──────────────────────────────────────────────────

greedy_forward <- function(pe, term_groups, min_events = 6L) {
  if (sum(pe$occurred) < min_events) return(NULL)

  null_fit <- tryCatch(
    coxph(Surv(time1, occurred) ~ strata(time2), data = pe),
    error = function(e) NULL
  )
  if (is.null(null_fit)) return(NULL)

  current_aic  <- AIC(null_fit)
  current_cols <- character(0)
  remaining    <- term_groups
  path         <- vector("list", length(term_groups))

  for (step in seq_along(term_groups)) {
    candidates <- vapply(names(remaining), function(gname) {
      cols <- c(current_cols, remaining[[gname]])
      # Drop columns absent from pe (guards against missing stats)
      cols <- cols[cols %in% names(pe)]
      if (length(cols) == 0) return(NA_real_)
      f <- as.formula(paste0(
        "Surv(time1, occurred) ~ ",
        paste(cols, collapse = " + "),
        " + strata(time2)"
      ))
      m <- tryCatch(coxph(f, data = pe), error = function(e) NULL)
      if (is.null(m)) NA_real_ else AIC(m)
    }, numeric(1))

    if (all(is.na(candidates))) break

    best_gname <- names(which.min(candidates))
    best_aic   <- candidates[best_gname]

    path[[step]] <- list(
      step           = step,
      added          = best_gname,
      aic            = best_aic,
      delta_aic      = best_aic - current_aic,
      prop_delta_aic = (best_aic - current_aic) / abs(current_aic),
      candidates     = candidates
    )

    current_cols <- c(current_cols, remaining[[best_gname]])
    remaining    <- remaining[names(remaining) != best_gname]
    current_aic  <- best_aic
  }

  path[!sapply(path, is.null)]
}

# ── Per-network fitting function ──────────────────────────────────────────────

fit_one_sep <- function(path, cat_label, visible_label) {
  e <- new.env(parent = emptyenv())
  tryCatch(load(path, envir = e), error = function(err) return(NULL))

  el_raw <- as.data.frame(get("edgelist", envir = e))[, 1:3]
  el     <- prepare_edgelist(el_raw)
  if (nrow(el) < 20L) return(NULL)

  all_nodes <- sort(unique(c(el$sender, el$receiver)))
  base_rs   <- llrem:::get_base_riskset(NULL, all_nodes)

  pe <- tryCatch(
    llrem:::remstats_generate(el, n = length(all_nodes),
                              n_nonevents = 200L, base_riskset = base_rs),
    error = function(err) NULL
  )
  if (is.null(pe)) return(NULL)
  pe <- llrem:::add_survival_cols(pe)
  if (!"SndInitiation" %in% names(pe))
    pe$SndInitiation <- pe$RecInitiation

  pe_init <- pe[pe$SndInitiation == 0, ]
  pe_cont <- pe[pe$SndInitiation == 1, ]

  path_init <- greedy_forward(pe_init, term_groups_init)
  path_cont <- greedy_forward(pe_cont, term_groups_cont)

  if (is.null(path_init) && is.null(path_cont)) return(NULL)

  list(
    file       = path,
    cat        = cat_label,
    visible    = visible_label,
    n          = length(all_nodes),
    n_events   = nrow(el),
    n_init     = sum(pe_init$occurred),
    n_cont     = sum(pe_cont$occurred),
    path_init  = path_init,
    path_cont  = path_cont
  )
}

# ── Parallel execution ────────────────────────────────────────────────────────

n_cores <- min(parallel::detectCores() - 1L, 4L)
cl <- parallel::makeCluster(n_cores)
parallel::clusterCall(cl, function() {
  devtools::load_all("llrem", quiet = TRUE)
  library(survival)
})
parallel::clusterExport(cl, c("term_groups_init", "term_groups_cont",
                               "greedy_forward", "fit_one_sep", "DATA_DIR"))

cat(sprintf("Running separable greedy selection on %d networks using %d cores ...\n",
            nrow(all_files), n_cores))

set.seed(8214L)
raw_results <- parallel::clusterMap(
  cl, fit_one_sep,
  path          = all_files$path,
  cat_label     = all_files$cat,
  visible_label = all_files$visible,
  SIMPLIFY      = FALSE
)
parallel::stopCluster(cl)

greedy_sep <- raw_results[!sapply(raw_results, is.null)]
cat(sprintf("%d networks retained (dropped %d)\n\n",
            length(greedy_sep), nrow(all_files) - length(greedy_sep)))

save(greedy_sep, file = file.path(OUT_DIR, "greedy_sep_results.Rdata"))

# ── Summary helpers ───────────────────────────────────────────────────────────

cats    <- sapply(greedy_sep, `[[`, "cat")
vis     <- sapply(greedy_sep, `[[`, "visible")
vis_bin <- ifelse(vis == "Y", "High", "Low")
n_nets  <- length(greedy_sep)

# Build order/prop matrices for one regime
make_matrices <- function(path_key, term_groups) {
  gnames  <- names(term_groups)
  n_steps <- length(gnames)

  order_mat    <- matrix(NA_character_, nrow = n_nets, ncol = n_steps)
  prop_mat     <- matrix(NA_real_,      nrow = n_nets, ncol = n_steps)

  for (i in seq_len(n_nets)) {
    path <- greedy_sep[[i]][[path_key]]
    if (is.null(path)) next
    for (s in seq_along(path)) {
      if (s > n_steps) break
      order_mat[i, s] <- path[[s]]$added
      prop_mat[i, s]  <- path[[s]]$prop_delta_aic
    }
  }
  list(order = order_mat, prop = prop_mat, gnames = gnames, n_steps = n_steps)
}

make_freq_mat <- function(order_mat, gnames, n_steps, idx = seq_len(nrow(order_mat))) {
  fm <- matrix(0, nrow = length(gnames), ncol = n_steps,
               dimnames = list(gnames, paste0("Step", seq_len(n_steps))))
  for (i in idx)
    for (s in seq_len(n_steps)) {
      g <- order_mat[i, s]
      if (!is.na(g) && g %in% gnames)
        fm[g, s] <- fm[g, s] + 1
    }
  fm / length(idx)
}

modal_order <- function(order_mat, idx = seq_len(nrow(order_mat))) {
  apply(order_mat[idx, , drop = FALSE], 2, function(col) {
    tb <- sort(table(col[!is.na(col)]), decreasing = TRUE)
    if (length(tb) == 0) return(NA_character_)
    paste0(names(tb)[1], " (", round(tb[1] / sum(tb) * 100), "%)")
  })
}

print_regime <- function(label, path_key, term_groups) {
  cat(rep("=", 65), "\n", sep = "")
  cat(sprintf("%s\n", label))
  cat(rep("=", 65), "\n\n", sep = "")

  mx   <- make_matrices(path_key, term_groups)
  gnames  <- mx$gnames
  n_steps <- mx$n_steps
  order_mat <- mx$order
  prop_mat  <- mx$prop
  cum_mat   <- t(apply(prop_mat, 1, function(r) cumsum(ifelse(is.na(r), 0, r))))

  # Selection frequency
  freq_mat <- make_freq_mat(order_mat, gnames, n_steps)
  cat("--- Selection frequency by step (proportion of networks) ---\n")
  print(round(freq_mat, 2))
  cat("\n")

  # Mean proportional delta AIC
  mean_step <- colMeans(prop_mat, na.rm = TRUE)
  sd_step   <- apply(prop_mat, 2, sd, na.rm = TRUE)
  mean_cum  <- colMeans(cum_mat, na.rm = TRUE)
  cat("--- Mean proportional delta AIC by step ---\n")
  n_obs <- colSums(!is.na(prop_mat))
  if (any(n_obs > 0)) {
    print(data.frame(
      Step       = seq_len(n_steps),
      PerStep    = round(mean_step, 4),
      SD         = round(sd_step,   4),
      Cumulative = round(mean_cum,  4),
      N          = n_obs
    ), row.names = FALSE)
  } else {
    cat("  (no data)\n")
  }
  cat("\n")

  # Modal order by domain
  cat("--- Most common greedy selection order by domain ---\n")
  for (g in sort(unique(cats))) {
    idx <- which(cats == g)
    mo  <- modal_order(order_mat, idx)
    cat(sprintf("  %-18s (n=%3d):  %s\n", g, length(idx),
                paste(mo[!is.na(mo)], collapse = " -> ")))
  }
  cat("\n")

  # Return for plotting
  invisible(list(freq_mat = freq_mat, order_mat = order_mat,
                 prop_mat = prop_mat, cum_mat = cum_mat,
                 gnames = gnames, n_steps = n_steps))
}

# ── Tables ────────────────────────────────────────────────────────────────────

init_mx <- print_regime(
  "INITIATION MODEL — greedy selection (SndInitiation == 0)",
  "path_init", term_groups_init
)

cont_mx <- print_regime(
  "CONTINUATION MODEL — greedy selection (SndInitiation == 1)",
  "path_cont", term_groups_cont
)

# ── Head-to-head: step-1 selection by domain ─────────────────────────────────

cat(rep("=", 65), "\n", sep = "")
cat("FIRST SELECTED TERM BY DOMAIN: INITIATION vs CONTINUATION\n")
cat(rep("=", 65), "\n\n", sep = "")

for (g in sort(unique(cats))) {
  idx <- which(cats == g)
  mo_i <- modal_order(init_mx$order_mat, idx)[1]
  mo_c <- modal_order(cont_mx$order_mat, idx)[1]
  cat(sprintf("  %-18s  init: %-30s  cont: %s\n", g,
              ifelse(is.na(mo_i), "---", mo_i),
              ifelse(is.na(mo_c), "---", mo_c)))
}
cat("\n")

# ── Plots ─────────────────────────────────────────────────────────────────────

heatmap_theme <- theme_minimal(base_size = 10) +
  theme(panel.grid    = element_blank(),
        strip.text    = element_text(face = "bold"),
        axis.text.y   = element_text(size = 8),
        legend.position = "right")

make_long_df <- function(mx, regime_label, idx = seq_len(n_nets)) {
  fm <- make_freq_mat(mx$order_mat, mx$gnames, mx$n_steps, idx)
  df <- as.data.frame(fm)
  df$term <- rownames(fm)
  pivot_longer(df, -term, names_to = "Step", values_to = "freq") |>
    mutate(
      step_num = as.integer(sub("Step", "", Step)),
      regime   = regime_label,
      label    = ifelse(freq >= 0.05, sprintf("%.0f", freq * 100), "")
    )
}

# Overall heatmap: initiation vs continuation side by side
all_mss_i <- rowSums(init_mx$freq_mat *
  matrix(seq_len(init_mx$n_steps), nrow = init_mx$n_steps,
         ncol = init_mx$n_steps, byrow = TRUE))
all_mss_c <- rowSums(cont_mx$freq_mat *
  matrix(seq_len(cont_mx$n_steps), nrow = cont_mx$n_steps,
         ncol = cont_mx$n_steps, byrow = TRUE))

term_ord_i <- names(sort(all_mss_i))
term_ord_c <- names(sort(all_mss_c))

df_init <- make_long_df(init_mx, "Initiation") |>
  mutate(term = factor(term, levels = term_ord_i))
df_cont <- make_long_df(cont_mx, "Continuation") |>
  mutate(term = factor(term, levels = term_ord_c))

plot_heatmap <- function(df, title) {
  ggplot(df, aes(x = step_num, y = term, fill = freq)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = label), size = 2.6, colour = "white", fontface = "bold") +
    scale_fill_gradient(low = "white", high = "#2171b5", limits = c(0, 1),
                        name = "Prop.\nnetworks") +
    scale_x_continuous(breaks = seq_len(max(df$step_num)), expand = c(0, 0)) +
    labs(x = "Greedy step", y = NULL, title = title) +
    heatmap_theme
}

p_init <- plot_heatmap(df_init, "Initiation regime: term selection frequency")
p_cont <- plot_heatmap(df_cont, "Continuation regime: term selection frequency")

ggsave(file.path(OUT_DIR, "greedy_sep_init_heatmap.pdf"),
       p_init, width = 9, height = 5)
ggsave(file.path(OUT_DIR, "greedy_sep_cont_heatmap.pdf"),
       p_cont, width = 9, height = 5)

# Visibility split: initiation
idx_hi <- which(vis_bin == "High")
idx_lo <- which(vis_bin == "Low")

df_vis_init <- bind_rows(
  make_long_df(init_mx, sprintf("High-vis (n=%d)", length(idx_hi)), idx_hi),
  make_long_df(init_mx, sprintf("Low-vis  (n=%d)", length(idx_lo)), idx_lo)
) |> mutate(term = factor(term, levels = term_ord_i))

df_vis_cont <- bind_rows(
  make_long_df(cont_mx, sprintf("High-vis (n=%d)", length(idx_hi)), idx_hi),
  make_long_df(cont_mx, sprintf("Low-vis  (n=%d)", length(idx_lo)), idx_lo)
) |> mutate(term = factor(term, levels = term_ord_c))

make_facet_heatmap <- function(df, title) {
  ggplot(df, aes(x = step_num, y = term, fill = freq)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = label), size = 2.4, colour = "white", fontface = "bold") +
    facet_wrap(~regime, ncol = 2) +
    scale_fill_gradient(low = "white", high = "#2171b5", limits = c(0, 1),
                        name = "Prop.\nnetworks") +
    scale_x_continuous(breaks = seq_len(max(df$step_num)), expand = c(0, 0)) +
    labs(x = "Greedy step", y = NULL, title = title) +
    heatmap_theme
}

ggsave(file.path(OUT_DIR, "greedy_sep_vis_init.pdf"),
       make_facet_heatmap(df_vis_init, "Initiation: selection frequency by visibility"),
       width = 12, height = 5)

ggsave(file.path(OUT_DIR, "greedy_sep_vis_cont.pdf"),
       make_facet_heatmap(df_vis_cont, "Continuation: selection frequency by visibility"),
       width = 12, height = 5)

cat("Results saved to greedy_sep_results.Rdata\n")
cat("Plots: greedy_sep_init_heatmap.pdf, greedy_sep_cont_heatmap.pdf\n")
cat("       greedy_sep_vis_init.pdf, greedy_sep_vis_cont.pdf\n")
cat("\nDone.\n")
