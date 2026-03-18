# This test verifies mathematical inversion consistency:
# N obtained from paired_ss_opt_bin should return target power
# when evaluated with paired_power_opt_bin



test_that("paired_ss_opt_bin and paired_power_opt_bin are consistent", {

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

  expect_equal(pow$overall_power, 0.80, tolerance = 1e-3)
})


