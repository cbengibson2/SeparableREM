make_two_models <- function() {
  m1 <- fit_rem(test_edgelist, test_formula,
                n_nonevents = NULL, base_riskset = test_riskset)
  res <- run_llm_rem(
    test_edgelist, test_formula,
    query_fn     = make_query_roleplay_sender(mock_cfg),
    n_nonevents  = NULL, base_riskset = test_riskset, seed = 1L
  )
  list("Empirical" = m1, "LLM" = res$rem)
}

test_that("compare_rem_models returns a ggplot invisibly", {
  models <- make_two_models()
  p <- suppressWarnings(compare_rem_models(models, title = "Test"))
  expect_s3_class(p, "gg")
})

test_that("compare_rem_models data has one row per term per model", {
  models   <- make_two_models()
  p        <- suppressWarnings(compare_rem_models(models))
  n_terms  <- length(names(stats::coef(models[[1L]])))
  expect_equal(nrow(p$data), n_terms * length(models))
})

test_that("term_order sets factor levels correctly", {
  models <- make_two_models()
  ord    <- c("NTriadic", "IOviolation", "Recencylog",
              "RecInitiation", "Repetition", "Imitation")
  p      <- suppressWarnings(compare_rem_models(models, term_order = ord))
  expect_equal(levels(p$data$term), ord)
})

test_that("custom palette is applied to rendered point colours", {
  models  <- make_two_models()
  palette <- c("Empirical" = "#FF0000", "LLM" = "#0000FF")
  p       <- suppressWarnings(compare_rem_models(models, palette = palette))
  built   <- suppressWarnings(ggplot2::ggplot_build(p))
  colours <- toupper(unique(built$data[[3L]]$colour))
  expect_true(all(c("#FF0000", "#0000FF") %in% colours))
})

test_that("compare_rem_models errors on an unnamed model list", {
  models <- make_two_models()
  expect_snapshot(compare_rem_models(unname(models)), error = TRUE)
})
