# Extracted from test-power_bin_consistency.R:15

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "DxPower", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
res <- ss_opt_bin(
    se.e =  0.85, sp.e = 0.85,
    se.min =  0.70, sp.min =  0.65,
    alpha = 0.025,
    prev = 0.5,
    power = 0.80
  )
