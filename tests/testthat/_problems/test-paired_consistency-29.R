# Extracted from test-paired_consistency.R:29

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "DxPower", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
res <- paired_ss_opt_bin(
    alpha = 0.05,
    power = 0.80,
    se.c = 0.81, se.e = 0.90,
    sp.c = 0.66, sp.e = 0.80,
    psi.d = 0.09,
    psi.nd = 0.14,
    prev = 0.47
  )
N_star <- res$N_tot
pow <- paired_power_opt_bin(
    N = N_star,
    alpha = 0.05,
    se.c = 0.81, se.e = 0.90,
    sp.c = 0.66, sp.e = 0.80,
    psi.d = 0.09,
    psi.nd = 0.14,
    prev = 0.47
  )
