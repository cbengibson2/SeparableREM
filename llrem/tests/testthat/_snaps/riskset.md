# prepare_edgelist errors on fewer than 3 columns

    Code
      prepare_edgelist(data.frame(a = 1:3, b = 1:3))
    Condition
      Error:
      ! `raw` must have at least 3 columns: time, sender, receiver.

# fit_rem errors on a 0-row edgelist

    Code
      fit_rem(test_edgelist[0L, ], test_formula, base_riskset = test_riskset)
    Condition
      Error in `$<-.data.frame`:
      ! replacement has 1 row, data has 0

