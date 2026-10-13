# Shared synthetic data and config used across test files.

make_edgelist <- function(n_nodes = 6L, n_events = 20L, seed = 1776L) {
  set.seed(seed)
  out <- data.frame(time = integer(n_events), sender = integer(n_events),
                    receiver = integer(n_events))
  for (i in seq_len(n_events)) {
    s <- sample.int(n_nodes, 1L)
    r <- sample(setdiff(seq_len(n_nodes), s), 1L)
    out[i, ] <- list(i, s, r)
  }
  prepare_edgelist(out)
}

test_edgelist  <- make_edgelist()
test_nodes     <- get_nodes(test_edgelist)
test_riskset   <- get_base_riskset(NULL, test_nodes)
test_formula   <- survival::Surv(time = time1, event = occurred) ~
  Imitation + Repetition + RecInitiation + Recencylog + IOviolation + NTriadic +
  survival::strata(time2)

mock_cfg <- rem_cfg(
  provider       = llm_provider_mock(mock_strategy_highest_indegree),
  context        = "A small test network.",
  action_verb    = "contacted",
  actor_noun     = "Node",
  history_window = NULL
)
