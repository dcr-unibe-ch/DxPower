# This test verifies mathematical inversion consistency:
# N obtained from power_opt should return target power
# when evaluated with power_diag.



test_that("power_opt and power_diag are consistent", {

  res <- power_opt(
    se1 = 0.85, sp1 = 0.85,
    se0 = 0.70, sp0 = 0.65,
    alpha = 0.025,
    prev = 0.5,
    pwr = 0.80
  )

  N_star <- res$N_tot

  pow <- power_diag(
    N = N_star,
    se1 = 0.85,
    se0 = 0.70,
    sp1 = 0.85,
    sp0 = 0.65,
    alpha = 0.025,
    prev = 0.5
  )

  expect_equal(pow$overall, 0.80, tolerance = 1e-3)
})
