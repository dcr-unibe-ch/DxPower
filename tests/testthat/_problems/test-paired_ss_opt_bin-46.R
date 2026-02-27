# Extracted from test-paired_ss_opt_bin.R:46

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "DxPower", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
res <- paired_ss_opt_bin(
    alpha = 0.05,
    power = 0.80,
    se.c = 0.80, se.e = 0.90,
    sp.c = 0.70, sp.e = 0.75,
    psi.d = 0.10,
    psi.nd = 0.20,
    prev = 0.5,
    delta.sp = 0.05,
    type.sp = "ni"
  )
expect_s3_class(res, "dx_bin")
expect_true(res$N_tot > 0)
expect_equal(res$target_power, 0.80, tolerance = 0.02)
