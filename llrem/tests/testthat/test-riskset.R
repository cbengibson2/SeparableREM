test_that("prepare_edgelist returns a data frame with correct columns", {
  raw <- data.frame(t = 1:5, s = c(1,2,3,1,2), r = c(2,3,1,3,1))
  out <- prepare_edgelist(raw)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("time", "sender", "receiver", "eveid"), ignore.order = FALSE)
})

test_that("prepare_edgelist coerces columns to integer", {
  raw <- data.frame(t = 1:3, s = c(1,2,3), r = c(2,3,1))
  out <- prepare_edgelist(raw)
  expect_type(out$time,     "integer")
  expect_type(out$sender,   "integer")
  expect_type(out$receiver, "integer")
})

test_that("prepare_edgelist forms eveid correctly", {
  raw <- data.frame(t = 1:2, s = c(1L, 3L), r = c(2L, 1L))
  out <- prepare_edgelist(raw)
  expect_equal(out$eveid, c("1_2", "3_1"))
})

test_that("prepare_edgelist accepts matrix input", {
  mat <- matrix(c(1,1,2, 2,2,3), nrow = 2, byrow = TRUE)
  out <- prepare_edgelist(mat)
  expect_s3_class(out, "data.frame")
  expect_true("eveid" %in% names(out))
})

test_that("prepare_edgelist errors on fewer than 3 columns", {
  expect_snapshot(prepare_edgelist(data.frame(a = 1:3, b = 1:3)), error = TRUE)
})

test_that("fit_rem returns a coxph object", {
  m <- fit_rem(test_edgelist, test_formula,
               n_nonevents = NULL, base_riskset = test_riskset)
  expect_s3_class(m, "coxph")
})

test_that("fit_rem produces expected coefficient names", {
  m <- fit_rem(test_edgelist, test_formula,
               n_nonevents = NULL, base_riskset = test_riskset)
  expected <- c("Imitation", "Repetition", "RecInitiation",
                "Recencylog", "IOviolation", "NTriadic")
  # coef() includes strata terms; check that the core covariates are present
  expect_true(all(expected %in% names(stats::coef(m))))
})

test_that("fit_rem coefficients are all finite", {
  m <- fit_rem(test_edgelist, test_formula,
               n_nonevents = NULL, base_riskset = test_riskset)
  expect_true(all(is.finite(stats::coef(m))))
})

test_that("fit_rem errors on a 0-row edgelist", {
  expect_snapshot(
    fit_rem(test_edgelist[0L, ], test_formula, base_riskset = test_riskset),
    error = TRUE
  )
})

# Triadic closure (C++ implementation) ----------------------------------------

test_that("triadic_counts returns correct OSP/ISP/ITP/OTP for a known graph", {
  # Triangle: 1->2, 2->3, 1->3
  # Node 1 out-neighbours: {2, 3}
  # Node 2 out-neighbours: {3}
  # Node 3 out-neighbours: {}
  # Node 1 in-neighbours:  {}
  # Node 2 in-neighbours:  {1}
  # Node 3 in-neighbours:  {1, 2}
  snd <- c(1L, 2L, 1L)
  rec <- c(2L, 3L, 3L)
  adj <- llrem:::build_adjacency(snd, rec)

  # For dyad 1->2:
  #   OSP: intersect(out[1], out[2]) = intersect({2,3},{3}) = {3} -> 1
  #   ISP: intersect(in[1],  in[2])  = intersect({},  {1}) = {}  -> 0
  #   ITP: intersect(in[1],  out[2]) = intersect({},  {3}) = {}  -> 0
  #   OTP: intersect(out[1], in[2])  = intersect({2,3},{1}) = {} -> 0
  counts <- llrem:::triadic_counts(adj, 1L, 2L)
  expect_equal(as.integer(counts[, "OSP"]), 1L)
  expect_equal(as.integer(counts[, "ISP"]), 0L)
  expect_equal(as.integer(counts[, "ITP"]), 0L)
  expect_equal(as.integer(counts[, "OTP"]), 0L)

  # For dyad 2->3:
  #   OSP: intersect(out[2], out[3]) = intersect({3},{}) = {} -> 0
  #   ISP: intersect(in[2],  in[3])  = intersect({1},{1,2}) = {1} -> 1
  #   ITP: intersect(in[2],  out[3]) = intersect({1},{})  = {} -> 0
  #   OTP: intersect(out[2], in[3])  = intersect({3},{1,2}) = {} -> 0
  counts2 <- llrem:::triadic_counts(adj, 2L, 3L)
  expect_equal(as.integer(counts2[, "ISP"]), 1L)
  expect_equal(as.integer(counts2[, "OSP"]), 0L)
})

test_that("NTriadic from C++ matches igraph reference on small graph", {
  # Build a small edgelist with enough events for triadic stats
  el <- prepare_edgelist(data.frame(
    t = 1:6,
    s = c(1L, 2L, 3L, 1L, 2L, 4L),
    r = c(2L, 3L, 1L, 3L, 1L, 1L)
  ))
  rs <- llrem:::get_base_riskset(NULL, llrem:::get_nodes(el))

  # C++ path (current implementation)
  cpp_rs <- llrem:::riskset_cov2(el[1:5, ], el[6, ], n_nonevents = NULL,
                                  base_riskset = rs)

  # Reference: compute triadic closure with igraph directly
  past <- el[1:5, ]
  g    <- igraph::graph_from_data_frame(past[, c("sender", "receiver")],
                                         directed = TRUE)
  ref_osp <- integer(nrow(rs))
  ref_isp <- integer(nrow(rs))
  ref_itp <- integer(nrow(rs))
  ref_otp <- integer(nrow(rs))
  in_g <- which(rs$sender %in% igraph::V(g) & rs$receiver %in% igraph::V(g))
  lookup    <- unique(c(rs$sender[in_g], rs$receiver[in_g]))
  outneighs <- lapply(lookup, \(a) igraph::neighbors(g, a, mode = "out"))
  inneighs  <- lapply(lookup, \(a) igraph::neighbors(g, a, mode = "in"))
  names(outneighs) <- names(inneighs) <- as.character(lookup)
  for (i in in_g) {
    a <- as.character(rs$sender[i]);  b <- as.character(rs$receiver[i])
    ref_osp[i] <- length(intersect(outneighs[[a]], outneighs[[b]]))
    ref_isp[i] <- length(intersect(inneighs[[a]],  inneighs[[b]]))
    ref_itp[i] <- length(intersect(inneighs[[a]],  outneighs[[b]]))
    ref_otp[i] <- length(intersect(outneighs[[a]], inneighs[[b]]))
  }

  expect_equal(cpp_rs$OSP, ref_osp)
  expect_equal(cpp_rs$ISP, ref_isp)
  expect_equal(cpp_rs$ITP, ref_itp)
  expect_equal(cpp_rs$OTP, ref_otp)
})

# Disaggregated normalised triadic terms --------------------------------------

test_that("NOTPSnd / NITPSnd / NOSPSnd / NISPSnd are present and in [0, 1]", {
  el <- prepare_edgelist(data.frame(
    t = 1:8,
    s = c(1L, 2L, 3L, 1L, 2L, 3L, 4L, 1L),
    r = c(2L, 3L, 1L, 3L, 4L, 2L, 1L, 4L)
  ))
  rs  <- llrem:::get_base_riskset(NULL, llrem:::get_nodes(el))
  out <- llrem:::riskset_cov2(el[1:7, ], el[8, ], n_nonevents = NULL,
                               base_riskset = rs)
  for (col in c("NOTPSnd", "NITPSnd", "NOSPSnd", "NISPSnd")) {
    expect_true(col %in% names(out), label = paste(col, "present"))
    expect_true(all(out[[col]] >= 0 & out[[col]] <= 1),
                label = paste(col, "in [0,1]"))
  }
})

test_that("NOTPSnd is zero wherever OTP is zero", {
  el <- prepare_edgelist(data.frame(
    t = 1:8,
    s = c(1L, 2L, 3L, 1L, 2L, 3L, 4L, 1L),
    r = c(2L, 3L, 1L, 3L, 4L, 2L, 1L, 4L)
  ))
  rs  <- llrem:::get_base_riskset(NULL, llrem:::get_nodes(el))
  out <- llrem:::riskset_cov2(el[1:7, ], el[8, ], n_nonevents = NULL,
                               base_riskset = rs)
  expect_true(all(out$NOTPSnd[out$OTP == 0] == 0))
  expect_true(all(out$NOSPSnd[out$OSP == 0] == 0))
})

test_that("per-sender max normalisation: max NOTPSnd per sender is 0 or 1", {
  el <- prepare_edgelist(data.frame(
    t = 1:8,
    s = c(1L, 2L, 3L, 1L, 2L, 3L, 4L, 1L),
    r = c(2L, 3L, 1L, 3L, 4L, 2L, 1L, 4L)
  ))
  rs  <- llrem:::get_base_riskset(NULL, llrem:::get_nodes(el))
  out <- llrem:::riskset_cov2(el[1:7, ], el[8, ], n_nonevents = NULL,
                               base_riskset = rs)
  sender_maxes <- tapply(out$NOTPSnd, out$sender, max)
  expect_true(all(sender_maxes == 0 | sender_maxes == 1))
})
