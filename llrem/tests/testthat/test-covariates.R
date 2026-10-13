# Tests for rem_covariates and covariate merging

# ---- constructors and validation --------------------------------------------

test_that("cov_node_static constructs and validates", {
  d <- data.frame(node = 1:3, value = c(1, 2, 3))
  cv <- cov_node_static("x", d)
  expect_s3_class(cv, "cov_node_static")
  expect_s3_class(cv, "rem_cov")
  expect_equal(cv$name, "x")
  expect_equal(cv$role, "both")
})

test_that("cov_node_static errors on bad data", {
  expect_error(cov_node_static("x", data.frame(a = 1, b = 2)), "node")
})

test_that("cov_node_temporal constructs and pre-sorts by node + time_start", {
  d <- data.frame(node = c(1L, 1L, 2L), time_start = c(5, 1, 1),
                  value = c(10, 5, 20))
  cv <- cov_node_temporal("y", d)
  expect_s3_class(cv, "cov_node_temporal")
  expect_equal(cv$by_node[["1"]]$time_start, c(1, 5))
})

test_that("cov_dyad_static constructs and pre-computes eid key", {
  d <- data.frame(sender = c(1L, 2L), receiver = c(2L, 1L), value = c(0.5, 0.3))
  cv <- cov_dyad_static("dist", d)
  expect_s3_class(cv, "cov_dyad_static")
  expect_true(".eid" %in% names(cv$data))
  expect_equal(cv$data$.eid, c("1_2", "2_1"))
})

test_that("cov_dyad_temporal constructs and pre-splits by dyad", {
  d <- data.frame(sender = c(1L, 1L), receiver = c(2L, 2L),
                  time_start = c(1, 10), value = c(0.1, 0.9))
  cv <- cov_dyad_temporal("trade", d)
  expect_s3_class(cv, "cov_dyad_temporal")
  expect_true("1_2" %in% names(cv$by_dyad))
})

test_that("rem_covariates validates name uniqueness", {
  d <- data.frame(node = 1:2, value = 1:2)
  expect_error(
    rem_covariates(cov_node_static("x", d), cov_node_static("x", d)),
    "Duplicate"
  )
})

test_that("rem_covariates blocks reserved column names", {
  d <- data.frame(node = 1:2, value = 1:2)
  expect_error(rem_covariates(cov_node_static("Repetition", d)), "reserved")
})

test_that("rem_covariates with role=both blocks <name>_snd collision", {
  d <- data.frame(node = 1:2, value = 1:2)
  # _snd suffix would collide if two covariates expand to same name
  cv1 <- cov_node_static("gdp", d, role = "sender")   # column: gdp
  cv2 <- cov_node_static("gdp", d, role = "receiver")  # column: gdp -> duplicate
  expect_error(rem_covariates(cv1, cv2), "Duplicate")
})

test_that("empty rem_covariates is valid", {
  rc <- rem_covariates()
  expect_s3_class(rc, "rem_covariates")
  expect_length(rc$covariates, 0L)
})

# ---- attach_covariates merge correctness ------------------------------------

make_rs <- function() {
  data.frame(
    sender   = c(1L, 1L, 2L, 2L),
    receiver = c(2L, 3L, 1L, 3L),
    eid      = c("1_2", "1_3", "2_1", "2_3"),
    stringsAsFactors = FALSE
  )
}

test_that("cov_node_static role=both adds _snd and _rec columns", {
  rs <- make_rs()
  d  <- data.frame(node = 1:3, value = c(10, 20, 30))
  rc <- rem_covariates(cov_node_static("pop", d, warn_missing = FALSE))
  out <- llrem:::attach_covariates(rs, rc, event_time = 1)
  expect_true("pop_snd" %in% names(out))
  expect_true("pop_rec" %in% names(out))
  # row 1: sender=1 -> pop_snd=10, receiver=2 -> pop_rec=20
  expect_equal(out$pop_snd, c(10, 10, 20, 20))
  expect_equal(out$pop_rec, c(20, 30, 10, 30))
})

test_that("cov_node_static role=sender adds one column named <name>", {
  rs  <- make_rs()
  d   <- data.frame(node = 1:3, value = c(1, 2, 3))
  rc  <- rem_covariates(cov_node_static("x", d, role = "sender",
                                         warn_missing = FALSE))
  out <- llrem:::attach_covariates(rs, rc, event_time = 1)
  expect_true("x" %in% names(out))
  expect_false("x_snd" %in% names(out))
  expect_equal(out$x, c(1, 1, 2, 2))
})

test_that("missing node match keeps NA and warns", {
  rs <- make_rs()
  d  <- data.frame(node = 1L, value = 5)   # only node 1 present
  rc <- rem_covariates(cov_node_static("x", d, role = "sender"))
  expect_warning(
    out <- llrem:::attach_covariates(rs, rc, event_time = 1),
    "no node match"
  )
  # node 2 has no match -> NA preserved
  expect_equal(out$x, c(5, 5, NA_real_, NA_real_))
})

test_that("cov_node_temporal step-function respects time windows", {
  rs <- make_rs()
  d  <- data.frame(
    node       = c(1L, 1L, 2L, 2L, 3L),
    time_start = c(1,  10,  1,  10,  1),
    value      = c(1,   2,  3,   4,  5)
  )
  rc <- rem_covariates(cov_node_temporal("x", d, role = "sender",
                                          warn_missing = FALSE))
  # event at t=5: nodes 1 and 2 in window 1 (time_start=1)
  out5 <- llrem:::attach_covariates(rs, rc, event_time = 5)
  expect_equal(out5$x, c(1, 1, 3, 3))

  # event at t=15: nodes 1 and 2 in window 2 (time_start=10)
  out15 <- llrem:::attach_covariates(rs, rc, event_time = 15)
  expect_equal(out15$x, c(2, 2, 4, 4))
})

test_that("cov_dyad_static fills matched values and keeps NA for unmatched", {
  rs <- make_rs()
  d  <- data.frame(sender = c(1L, 2L), receiver = c(2L, 3L), value = c(100, 200))
  rc <- rem_covariates(cov_dyad_static("dist", d, warn_missing = FALSE))
  out <- llrem:::attach_covariates(rs, rc, event_time = 1)
  # 1_2 -> 100, 1_3 -> NA, 2_1 -> NA, 2_3 -> 200
  expect_equal(out$dist, c(100, NA_real_, NA_real_, 200))
})

test_that("cov_dyad_temporal step-function lookup is correct", {
  rs <- make_rs()
  d  <- data.frame(
    sender     = c(1L, 1L),
    receiver   = c(2L, 2L),
    time_start = c(1,  10),
    value      = c(0.1, 0.9)
  )
  rc <- rem_covariates(cov_dyad_temporal("trade", d, warn_missing = FALSE))
  out7  <- llrem:::attach_covariates(rs, rc, event_time = 7)
  out12 <- llrem:::attach_covariates(rs, rc, event_time = 12)
  expect_equal(out7$trade[out7$eid == "1_2"],  0.1)
  expect_equal(out12$trade[out12$eid == "1_2"], 0.9)
})

# ---- integration: fit_rem with covariates -----------------------------------

test_that("fit_rem clogit with a node covariate runs and coefficient is finite", {
  el  <- make_edgelist(n_nodes = 6L, n_events = 25L, seed = 99L)
  rs  <- llrem:::get_base_riskset(NULL, llrem:::get_nodes(el))

  # simple static covariate: node "size" 1-6
  d   <- data.frame(node = 1:6, value = as.numeric(1:6))
  rc  <- rem_covariates(cov_node_static("size", d, warn_missing = FALSE))

  f   <- survival::Surv(time1, occurred) ~
    Repetition + Recencylog + IOviolation + size_snd + size_rec +
    survival::strata(time2)

  suppressWarnings(
    m <- fit_rem(el, f, n_nonevents = NULL, base_riskset = rs,
                 method = "clogit", covariates = rc)
  )
  expect_true(is.finite(m$coef["size_snd"]))
  expect_true(is.finite(m$coef["size_rec"]))
})

test_that("fit_rem clogit with a dyad covariate runs and coefficient is finite", {
  el    <- make_edgelist(n_nodes = 6L, n_events = 25L, seed = 42L)
  rs    <- llrem:::get_base_riskset(NULL, llrem:::get_nodes(el))
  nodes <- llrem:::get_nodes(el)

  # dyad covariate: random distance between all pairs
  set.seed(7)
  dyad_d <- expand.grid(sender = nodes, receiver = nodes)
  dyad_d <- dyad_d[dyad_d$sender != dyad_d$receiver, ]
  dyad_d$value <- runif(nrow(dyad_d))
  rc  <- rem_covariates(cov_dyad_static("dist", dyad_d, warn_missing = FALSE))

  f   <- survival::Surv(time1, occurred) ~
    Repetition + dist + survival::strata(time2)

  suppressWarnings(
    m <- fit_rem(el, f, n_nonevents = NULL, base_riskset = rs,
                 method = "clogit", covariates = rc)
  )
  expect_true(is.finite(m$coef["dist"]))
})

test_that("endogenous NAs are filled to 0; exogenous NAs are preserved", {
  el  <- make_edgelist(n_nodes = 6L, n_events = 15L, seed = 5L)
  rs  <- llrem:::get_base_riskset(NULL, llrem:::get_nodes(el))
  # covariate only covers nodes 1-3, leaving nodes 4-6 as NA
  d   <- data.frame(node = 1:3, value = c(10, 20, 30))
  rc  <- rem_covariates(cov_node_static("pop", d, role = "sender",
                                         warn_missing = FALSE))
  suppressWarnings(
    pe <- llrem:::remstats_generate(el, n_nonevents = NULL,
                                    base_riskset = rs, covariates = rc)
  )
  # endogenous: Repetition should have no NAs
  expect_false(anyNA(pe$Repetition))
  # exogenous: pop should have NAs for senders 4-6
  expect_true(anyNA(pe$pop))
  expect_true(all(!is.na(pe$pop[pe$sender %in% 1:3])))
  expect_true(all( is.na(pe$pop[pe$sender %in% 4:6])))
})

test_that("NULL covariates produces identical output to no covariates argument", {
  el <- make_edgelist(n_nodes = 5L, n_events = 15L, seed = 1L)
  rs <- llrem:::get_base_riskset(NULL, llrem:::get_nodes(el))
  f  <- survival::Surv(time1, occurred) ~
    Repetition + Recencylog + survival::strata(time2)

  suppressWarnings({
    m1 <- fit_rem(el, f, n_nonevents = NULL, base_riskset = rs, method = "clogit")
    m2 <- fit_rem(el, f, n_nonevents = NULL, base_riskset = rs, method = "clogit",
                  covariates = NULL)
    m3 <- fit_rem(el, f, n_nonevents = NULL, base_riskset = rs, method = "clogit",
                  covariates = rem_covariates())
  })
  expect_equal(m1$coef, m2$coef)
  expect_equal(m1$coef, m3$coef)
})
