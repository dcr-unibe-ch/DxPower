# Extracted from test-power_opt_bin.R:11

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "DxPower", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
pow <- power_opt_bin(
    N=133,
    se.e = 0.85,
    se.min = 0.70,
    sp.e = 0.85,
    sp.min = 0.65,
    alpha = 0.025,
    prev = 0.5
  )
