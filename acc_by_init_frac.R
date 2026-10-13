# acc_by_init_frac.R
#
# Plots proportional improvement in predictive accuracy (separable vs pooled)
# against initiation fraction, one point per network.
# Uses acc_sep_cont vs acc_pool (continuation sub-model accuracy, available for
# all networks). Points are coloured by domain; domain means are overlaid.

rm(list = ls())

load("revision_results.Rdata")
library(ggplot2)

# Proportional improvement: (sep_cont - pool) / pool
res_df$init_frac    <- res_df$n_init / res_df$n_events
res_df$prop_improve <- (res_df$acc_sep_cont - res_df$acc_pool) / res_df$acc_pool

df <- res_df[is.finite(res_df$prop_improve) & is.finite(res_df$init_frac), ]

# Domain means
domain_means <- aggregate(
  cbind(init_frac, prop_improve) ~ cat,
  data = df, FUN = mean
)

# Consistent domain colours
domain_levels <- c("Glances","Email","War","SMS","Classroom",
                   "Twitter","Phone (DoD)","PMs","Gang Violence")
domain_colours <- setNames(
  c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
    "#EDC948","#B07AA1","#FF9DA7","#9C755F"),
  domain_levels
)

df$cat            <- factor(df$cat,            levels = domain_levels)
domain_means$cat  <- factor(domain_means$cat,  levels = domain_levels)

p <- ggplot(df, aes(x = init_frac, y = prop_improve, colour = cat)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_point(alpha = 0.35, size = 1.6) +
  geom_point(data = domain_means, aes(x = init_frac, y = prop_improve, colour = cat),
             shape = 18, size = 4.5) +
  geom_smooth(data = df, aes(x = init_frac, y = prop_improve),
              inherit.aes = FALSE,
              method = "loess", span = 0.9,
              colour = "grey30", fill = "grey80", alpha = 0.3,
              linewidth = 0.7) +
  geom_text(data = domain_means, aes(label = cat, colour = cat),
            size = 2.8, fontface = "bold", vjust = -0.9, show.legend = FALSE) +
  scale_colour_manual(values = domain_colours, guide = "none") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     name = "Initiation fraction") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     name = "Proportional improvement in predictive accuracy\n(separable cont. vs pooled)") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9, colour = "grey40")
  ) +
  labs(
    title = "Separable model predictive gain by initiation fraction",
    subtitle = "Small points: individual networks. Diamonds: domain means. Shaded band: loess."
  )

out_path <- "acc_by_init_frac.pdf"
ggsave(out_path, p, width = 7, height = 5)
cat(sprintf("Written: %s\n", out_path))

# Console summary
cat("\nDomain means:\n")
domain_means <- domain_means[order(domain_means$init_frac), ]
cat(sprintf("%-20s  init_frac=%5.1f%%  prop_improve=%+.1f%%\n",
    domain_means$cat,
    100 * domain_means$init_frac,
    100 * domain_means$prop_improve))
