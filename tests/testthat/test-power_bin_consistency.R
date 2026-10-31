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

  N_star <- res$n

  pow <- power_opt_bin(
    N = N_star,
    se.e = 0.85,
    se.min =  0.70,
    sp.e =0.85,
    sp.min = 0.65,
    alpha = 0.025,
    prev = 0.5
  )

  # n is the smallest integer reaching the target power
  expect_equal(N_star, round(N_star))
  expect_gte(pow$power, 0.80)
  expect_lt(power_opt_bin(
    N = N_star - 1,
    se.e = 0.85, se.min = 0.70,
    sp.e = 0.85, sp.min = 0.65,
    alpha = 0.025, prev = 0.5)$power, 0.80)
})

test_that("ss_opt_bin works at extreme prevalence and with logit/two-sided", {
  for (prev in c(0.1, 0.9)) {
    for (ci in c("wald", "logit")) {
      res <- ss_opt_bin(power = 0.8, prev = prev,
                        se.e = 0.85, se.min = 0.70,
                        sp.e = 0.85, sp.min = 0.65,
                        ci_method = ci, alternative = "two.sided")
      expect_false(is.na(res$n))
      expect_gte(res$power, 0.8)
    }
  }
})

test_that("overall power equals product of endpoint powers", {
  args <- list(N = 130, prev = 0.5, se.e = 0.85, se.min = 0.70,
               sp.e = 0.85, sp.min = 0.65)
  p0 <- do.call(power_opt_bin, args)
  expect_equal(p0$power, p0$power_se * p0$power_sp, tolerance = 1e-10)
})


