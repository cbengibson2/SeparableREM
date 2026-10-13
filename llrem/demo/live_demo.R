# demo/live_demo.R
#
# llrem pipeline with real LLM providers.
#
# Before running:
#   - Set the relevant API key environment variable (see Section 1).
#   - Swap the provider constructor and model name to match your setup.
#   - Adjust n_events and n_nonevents for cost/speed trade-offs.
#
# The script is written so you can run it section by section interactively.

library(llrem)
library(survival)

# ── 1. Choose a provider ──────────────────────────────────────────────────────
# Uncomment exactly one block, or adapt llm_provider_openai_compat() for any
# OpenAI-wire-format endpoint (Together AI, Fireworks, Azure OpenAI, etc.).

# OpenAI (default gpt-4o)
# Sys.setenv(OPENAI_API_KEY = "sk-...")
provider <- llm_provider_openai()

# Anthropic Claude
# Sys.setenv(ANTHROPIC_API_KEY = "sk-ant-...")
# provider <- llm_provider_anthropic(model = "claude-3-5-sonnet-20241022")

# Google Gemini
# Sys.setenv(GEMINI_API_KEY = "AIza...")
# provider <- llm_provider_gemini(model = "gemini-2.0-flash")

# xAI Grok
# Sys.setenv(XAI_API_KEY = "xai-...")
# provider <- llm_provider_grok(model = "grok-3")

# Local Ollama (no key required; Ollama must be running)
# provider <- llm_provider_ollama(model = "llama3.2")

# Generic OpenAI-compatible endpoint
# Sys.setenv(MY_API_KEY = "...")
# provider <- llm_provider_openai_compat(
#   api_key  = Sys.getenv("MY_API_KEY"),
#   model    = "my-model-name",
#   base_url = "https://my.endpoint/v1/chat/completions",
#   name     = "MyProvider"
# )

print(provider)

# ── 2. Empirical edgelist ─────────────────────────────────────────────────────
# Replace this block with your own data.
#
# prepare_edgelist() expects columns in order: time, sender, receiver.
# Node IDs must be positive integers. Smaller sequences run faster.
#
# Example with a CSV:
#   raw      <- read.csv("data/my_events.csv")   # must have time, sender, receiver cols
#   edgelist <- prepare_edgelist(raw[, c("time", "sender", "receiver")])
#
# Synthetic stand-in (replace before running a real study):

set.seed(6047)
n_nodes  <- 10L
n_events <- 80L   # trim for speed; 200-500 is a realistic study size

make_demo_edgelist <- function(n_nodes, n_events) {
  rows <- vector("list", n_events)
  for (i in seq_len(n_events)) {
    s <- sample.int(n_nodes, 1L)
    r <- sample(setdiff(seq_len(n_nodes), s), 1L)
    rows[[i]] <- data.frame(time = i, sender = s, receiver = r)
  }
  prepare_edgelist(do.call(rbind, rows))
}

edgelist <- make_demo_edgelist(n_nodes, n_events)
cat("Events:", nrow(edgelist), "| Nodes:", n_distinct_nodes <- length(unique(c(edgelist$sender, edgelist$receiver))), "\n")

# ── 3. Node labels (optional) ─────────────────────────────────────────────────
# If your nodes have meaningful names, supply them here so the LLM receives
# human-readable labels in its prompts (e.g., country names, person names).
# Leave as NULL to use generic "Node 1", "Node 2", … labels.

node_labels <- NULL
# Example:
# node_labels <- c("1" = "USA", "2" = "Russia", "3" = "China",
#                  "4" = "UK",  "5" = "France", "6" = "Germany")

# ── 4. Configuration ──────────────────────────────────────────────────────────
# Tailor context, action_verb, and actor_noun to your domain.

cfg <- rem_cfg(
  provider       = provider,
  context        = paste(
    "You are participating in a directed communication network.",
    "Nodes send messages to one another over time."
  ),
  action_verb    = "messaged",
  actor_noun     = "Node",
  node_labels    = node_labels,
  history_window = 30L   # number of past events shown in each prompt
)

print(cfg)

# ── 5. Shared objects ─────────────────────────────────────────────────────────
# Pre-compute base risk set once; reuse across all fit_rem() / run_llm_rem() calls.

all_nodes    <- sort(unique(c(edgelist$sender, edgelist$receiver)))
base_riskset <- llrem:::get_base_riskset(NULL, all_nodes)

rem_formula <- Surv(time = time1, event = occurred) ~
  Imitation + Repetition + RecInitiation + Recencylog + IOviolation + NTriadic +
  strata(time2)

# ── 6. Empirical REM ──────────────────────────────────────────────────────────

cat("\nFitting empirical REM...\n")
rem_empirical <- fit_rem(
  edgelist,
  formula      = rem_formula,
  n_nonevents  = 200L,
  base_riskset = base_riskset
)
print(summary(rem_empirical)$coefficients)

# ── 7. LLM sequence generation ────────────────────────────────────────────────
# Choose one or more query strategies. Each call hits the API n_events times,
# so budget accordingly (~$0.01–0.10 per event depending on model and context).
#
# Checkpointing writes progress every 10 events; remove the file to restart
# from scratch, or leave it to resume an interrupted run.

checkpoint_dir <- "demo/checkpoints"
dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)

# --- roleplay_sender: LLM adopts the empirical sender's identity ---
cat("\nRunning roleplay_sender...\n")
res_sender <- run_llm_rem(
  edgelist        = edgelist,
  formula         = rem_formula,
  query_fn        = make_query_roleplay_sender(cfg),
  n_nonevents     = 200L,
  base_riskset    = base_riskset,
  seed            = 9173L,
  checkpoint_path = file.path(checkpoint_dir, "sender_checkpoint.RData"),
  rate_limit_secs = 0.5   # be polite to the API; increase if you hit 429s
)

cat("Generated sequence (first 6 rows):\n")
print(head(res_sender$sequence))

# --- predict: LLM picks both sender and receiver ---
cat("\nRunning predict...\n")
res_predict <- run_llm_rem(
  edgelist        = edgelist,
  formula         = rem_formula,
  query_fn        = make_query_predict(cfg),
  n_nonevents     = 200L,
  base_riskset    = base_riskset,
  seed            = 1384L,
  checkpoint_path = file.path(checkpoint_dir, "predict_checkpoint.RData"),
  rate_limit_secs = 0.5
)

# ── 8. Inspect raw LLM responses ──────────────────────────────────────────────
# res_sender$responses is a data frame with one row per generated event, useful
# for auditing parse failures and understanding model behavior.

cat("\nRaw LLM response sample (sender run):\n")
print(head(res_sender$responses))

# ── 9. Forest-plot comparison ─────────────────────────────────────────────────

compare_rem_models(
  models = list(
    Empirical       = rem_empirical,
    `LLM (sender)`  = res_sender$rem,
    `LLM (predict)` = res_predict$rem
  ),
  title = paste("REM Coefficient Comparison —", provider$name, provider$model)
)

# ── 10. Save results ──────────────────────────────────────────────────────────

save(rem_empirical, res_sender, res_predict,
     file = file.path(checkpoint_dir, "results.RData"))
cat("\nResults saved to", file.path(checkpoint_dir, "results.RData"), "\n")
