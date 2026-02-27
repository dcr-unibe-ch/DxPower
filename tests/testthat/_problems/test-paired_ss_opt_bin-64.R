# Extracted from test-paired_ss_opt_bin.R:64

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "DxPower", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
ss_sup <- paired_ss_opt_bin(
    alpha = 0.05,
    power = 0.80,
    se.c = 0.80, se.e = 0.90,
    sp.c = 0.70, sp.e = 0.75,
    psi.d = 0.10,
    psi.nd = 0.20,
    prev = 0.5
  )
