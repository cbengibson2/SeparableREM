#' Forest plot comparing fitted REM coefficients
#'
#' Produces a coefficient plot with 95% confidence intervals for a named list
#' of `coxph` models, with one colour per model source.
#'
#' @param models Named list of `coxph` objects.
#' @param term_order Character vector giving the desired y-axis order (bottom
#'   to top). Defaults to reverse alphabetical over all terms present.
#' @param palette Named character vector of hex colours keyed to
#'   `names(models)`. A colour-blind-friendly default palette is applied when
#'   `NULL`.
#' @param title Plot title string.
#'
#' @return A `ggplot` object (invisibly). The plot is also printed as a
#'   side-effect.
#' @export
compare_rem_models <- function(models,
                                term_order = NULL,
                                palette    = NULL,
                                title      = "REM Coefficient Comparison") {
  if (!is.list(models) || is.null(names(models)))
    stop("`models` must be a named list of coxph objects.", call. = FALSE)

  extract_one <- function(model, label) {
    s <- summary(model)$coefficients
    data.frame(term    = rownames(s),
               coef    = s[, "coef"],
               se      = s[, "se(coef)"],
               p       = s[, "Pr(>|z|)"],
               lower95 = s[, "coef"] - 1.96 * s[, "se(coef)"],
               upper95 = s[, "coef"] + 1.96 * s[, "se(coef)"],
               source  = label,
               row.names = NULL, stringsAsFactors = FALSE)
  }

  coef_df <- do.call(rbind, Map(extract_one, models, names(models)))

  if (is.null(term_order)) term_order <- rev(sort(unique(coef_df$term)))
  if (is.null(palette)) {
    defaults <- c("#2166AC", "#D6604D", "#4DAC26", "#A858A0",
                  "#F4A736", "#878787", "#01665E", "#8C510A")
    palette  <- stats::setNames(defaults[seq_along(models)], names(models))
  }

  coef_df$term   <- factor(coef_df$term,   levels = term_order)
  coef_df$source <- factor(coef_df$source, levels = names(models))

  p <- ggplot2::ggplot(
    coef_df,
    ggplot2::aes(x = coef, y = term, colour = source,
                 xmin = lower95, xmax = upper95)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_linerange(
      position = ggplot2::position_dodge(width = 0.6), linewidth = 0.8
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.6), size = 2.5
    ) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::labs(title = title,
                  x = "Coefficient (log hazard ratio)", y = NULL, colour = NULL) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "bottom")

  print(p)
  invisible(p)
}
