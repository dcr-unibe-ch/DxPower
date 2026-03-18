# Extracted from test-paired_power_opt_bin.R:56

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "DxPower", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
p_small <- paired_power_opt_bin(
    N = 80,
    alpha = 0.05,
    se.c = 0.81, se.e = 0.90,
    sp.c = 0.66, sp.e = 0.80,
    psi.d = 0.09,
    psi.nd = 0.14,
    prev = 0.47
  )
