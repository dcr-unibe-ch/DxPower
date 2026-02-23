# This test verifies mathematical inversion consistency:
# N obtained from ss_opt_bin should return target power
# when evaluated with power_opt_bin



test_that("ss_opt_bin and power_opt_bin are consistent", {

  res <- ss_opt_bin(
    se.e =  0.85, sp.e = 0.85,
    se.min =  0.70, sp.min =  0.65,
    alpha = 0.025,
    prev = 0.5,
    power = 0.80
  )

  N_star <- res$N_tot

  pow <- power_opt_bin(
    N = N_star,
    se.e = 0.85,
    se.min =  0.70,
    sp.e =0.85,
    sp.min = 0.65,
    alpha = 0.025,
    prev = 0.5
  )

  expect_equal(pow$overall, 0.80, tolerance = 1e-3)
})
