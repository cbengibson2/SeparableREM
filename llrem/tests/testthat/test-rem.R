test_that("run_llm_rem (roleplay_sender) returns correct structure", {
  res <- run_llm_rem(
    test_edgelist, test_formula,
    query_fn     = make_query_roleplay_sender(mock_cfg),
    n_nonevents  = NULL, base_riskset = test_riskset, seed = 1L
  )
  expect_named(res, c("sequence", "responses", "rem"))
  expect_s3_class(res$rem, "coxph")
})

test_that("run_llm_rem sequence has the correct number of rows", {
  res <- run_llm_rem(
    test_edgelist, test_formula,
    query_fn     = make_query_roleplay_sender(mock_cfg),
    n_nonevents  = NULL, base_riskset = test_riskset, seed = 1L
  )
  expect_equal(nrow(res$sequence), nrow(test_edgelist))
})

test_that("roleplay_sender preserves the empirical sender sequence", {
  res <- run_llm_rem(
    test_edgelist, test_formula,
    query_fn     = make_query_roleplay_sender(mock_cfg),
    n_nonevents  = NULL, base_riskset = test_riskset, seed = 1L
  )
  expect_identical(res$sequence$sender, test_edgelist$sender)
})

test_that("generated sequence contains no self-loops", {
  res <- run_llm_rem(
    test_edgelist, test_formula,
    query_fn     = make_query_roleplay_sender(mock_cfg),
    n_nonevents  = NULL, base_riskset = test_riskset, seed = 1L
  )
  expect_true(all(res$sequence$sender != res$sequence$receiver))
})

test_that("all generated receivers are valid node IDs", {
  res <- run_llm_rem(
    test_edgelist, test_formula,
    query_fn     = make_query_roleplay_sender(mock_cfg),
    n_nonevents  = NULL, base_riskset = test_riskset, seed = 1L
  )
  expect_true(all(res$sequence$receiver %in% test_nodes))
})

test_that("run_llm_rem (predict) produces no self-loops and valid IDs", {
  res <- run_llm_rem(
    test_edgelist, test_formula,
    query_fn     = make_query_predict(mock_cfg),
    n_nonevents  = NULL, base_riskset = test_riskset, seed = 2L
  )
  expect_true(all(res$sequence$sender   %in% test_nodes))
  expect_true(all(res$sequence$receiver %in% test_nodes))
  expect_true(all(res$sequence$sender   != res$sequence$receiver))
})

test_that("run_llm_rem (roleplay_random) senders deviate from empirical", {
  res <- run_llm_rem(
    test_edgelist, test_formula,
    query_fn     = make_query_roleplay_random(mock_cfg),
    n_nonevents  = NULL, base_riskset = test_riskset, seed = 3L
  )
  # Exclude the seeded first event; with 6 nodes at least one should differ
  n_diff <- sum(res$sequence$sender[-1L] != test_edgelist$sender[-1L])
  expect_gte(n_diff, 1L)
})
