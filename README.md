# SepREM

Replication code for the separable relational event model (SepREM) paper.

Preprint: https://osf.io/3v5e2

## Overview

Standard relational event models (REMs) fit a single proportional-hazards model to all events in a sequence, treating first contact between two actors the same as repeated contact. SepREM partitions the event history by whether sender *i* has previously sent to receiver *j* (initiation vs. continuation), fits two independent sub-models on the disjoint subsets, and shows via AIC comparison that separable estimation outperforms the pooled model across 9 event-sequence domains (classroom interaction, email, gang violence, private messages, phone calls, SMS, Twitter, war/MID, face-to-face glances).

The `llrem/` subdirectory is an R package implementing the REM likelihood, risk-set statistics, and LLM-based sending strategies used across this and related projects.

## Run order

```
1. revision.R              -> revision_results.Rdata   (primary model fits)
2. greedy_opt_sep.R        -> greedy_sep_results.Rdata (greedy variable selection paths)
3. greedy_acc2.R           -> greedy_acc2_results.Rdata (predictive accuracy)
4. patch_acc.R             (patches NAs from parallel workers in revision_results.Rdata)
5. replicate.R             (console audit: prints every claim in the paper)
6. coef_table.R            -> paper/coef_table_sep.tex, paper/coef_table_pool.tex
7. sign_table.R            -> paper/sign_table.tex
8. acc_by_init_frac.R      -> acc_by_init_frac.pdf
```

## Scripts

| Script | Description |
|--------|-------------|
| `separable_all.R` | Separable and pooled REM on all domains; AIC comparison |
| `revision.R` | Revision analyses: interacted model, decay-window robustness, predictive accuracy, ecology-weighted summaries |
| `greedy_opt_sep.R` | Greedy forward AIC selection for the separable model across all networks |
| `greedy_acc2.R` | Predictive accuracy for greedy-selected separable and pooled models |
| `patch_acc.R` | Patches missing accuracy columns in `revision_results.Rdata` (sequential fallback for fork-based parallel worker failures) |
| `replicate.R` | Replication audit: loads stored results and prints every paper claim with the computed value |
| `coef_table.R` | Produces LaTeX coefficient tables for separable and pooled models |
| `sign_table.R` | Produces sign+significance table across all statistics and domains |
| `acc_by_init_frac.R` | Plots predictive accuracy gain vs. initiation fraction, one point per network |

## Sub-package: `llrem/`

An R package providing:
- `fit_rem()` — fits a REM via Cox proportional hazards or conditional logit
- `run_llm_rem()` — wraps an LLM as a sending-decision agent in a relational event simulation
- Risk-set statistics (participation shifts, recency, triadic terms, etc.)
- Demo scripts in `llrem/demo/`

Install locally with `devtools::install("llrem")` or load in-place with `devtools::load_all("llrem")`.

## Data

Network data (`DataP/`) is not included. Set the `DATA_DIR` environment variable to the directory containing domain subdirectories (`Classroom/`, `Email/`, `Gangs/`, etc.) before running any analysis script:

```r
Sys.setenv(DATA_DIR = "/path/to/DataP")
source("revision.R")
```

## Requirements

- R (≥ 4.1)
- Packages: `devtools`, `survival`, `parallel`, `ggplot2`, `tidyr`, `dplyr`
- A C++ compiler for `llrem`'s Rcpp components (`llrem/src/triadic.cpp`)

## Citation

See preprint at https://osf.io/3v5e2 for citation information.
