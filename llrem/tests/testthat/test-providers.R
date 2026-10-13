test_that("llm_provider_mock returns llm_provider class", {
  p <- llm_provider_mock(mock_strategy_min_id)
  expect_s3_class(p, "llm_provider")
})

test_that("print.llm_provider is silent and returns invisibly", {
  p <- llm_provider_mock(mock_strategy_min_id)
  expect_output(print(p), "Mock")
})

test_that("llm_call errors on non-provider input", {
  expect_snapshot(llm_call("not_a_provider", list()), error = TRUE)
})

test_that("mock_strategy_min_id returns parseable string", {
  p <- llm_provider_mock(mock_strategy_min_id)
  msg <- list(
    list(role = "system", content = "Test."),
    list(role = "user", content = paste0(
      "=== EVENT HISTORY ===\n",
      "  Event 1: Node 1 (1) contacted Node 3 (3).\n",
      "\n=== YOUR TURN ===\n",
      "You are Node 2 (ID: 2).\n",
      "Valid targets: Node 1 (1), Node 3 (3), Node 4 (4)\n",
      "Reply with the integer ID only."
    ))
  )
  resp <- p$call(msg)
  expect_type(resp, "character")
  expect_false(is.na(resp))
  expect_true(as.integer(resp) %in% c(1L, 3L, 4L))
})

test_that("mock_strategy_highest_indegree picks most-targeted node", {
  p <- llm_provider_mock(mock_strategy_highest_indegree)
  # Node 3 receives 2 events, node 1 receives 1 — expect "3"
  msg <- list(
    list(role = "system", content = "Test."),
    list(role = "user", content = paste0(
      "=== EVENT HISTORY ===\n",
      "  Event 1: Node 1 (1) contacted Node 3 (3).\n",
      "  Event 2: Node 2 (2) contacted Node 3 (3).\n",
      "  Event 3: Node 3 (3) contacted Node 1 (1).\n",
      "\n=== YOUR TURN ===\n",
      "You are Node 2 (ID: 2).\n",
      "Valid targets: Node 1 (1), Node 3 (3), Node 4 (4)\n",
      "Reply with the integer ID only."
    ))
  )
  expect_equal(p$call(msg), "3")
})

test_that("mock predict strategy returns two distinct valid integers", {
  p <- llm_provider_mock(mock_strategy_highest_indegree)
  msg <- list(
    list(role = "system", content = "Predict. Reply: SENDER_ID RECEIVER_ID"),
    list(role = "user", content = paste0(
      "=== EVENT HISTORY ===\n",
      "  Event 1: Node 1 (1) contacted Node 3 (3).\n",
      "\n=== PREDICTION TASK ===\n",
      "Actors: Node 1 (1), Node 2 (2), Node 3 (3), Node 4 (4)\n",
      "Reply: SENDER_ID RECEIVER_ID"
    ))
  )
  resp  <- p$call(msg, max_tokens = 25L)
  parts <- suppressWarnings(as.integer(strsplit(trimws(resp), "\\s+")[[1L]]))
  expect_length(parts, 2L)
  expect_false(anyNA(parts))
  expect_false(parts[1L] == parts[2L])
  expect_true(all(parts %in% 1:4))
})
