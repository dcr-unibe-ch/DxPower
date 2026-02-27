test_that("paired_power_opt_bin compute power for a fixed smaple size for sup test", {

  pow <- paired_power_opt_bin(
    N = 134,
    alpha = 0.05,
    se.c = 0.81, se.e = 0.90,
    sp.c = 0.66, sp.e = 0.80,
    psi.d = 0.09,
    psi.nd = 0.14,
    prev = 0.47
  )


  expect_s3_class(pow, "dx_bin")
  expect_true(pow$overall_power >= 0)
  expect_true(pow$overall_power <= 1)

  expect_true(all(c("power_se","power_sp","N_se","N_sp","N_tot") %in% names(pow)))
})


test_that("paired_power_opt_bin compute power for a fixed sample size for non-inferiority", {

  pow <- paired_power_opt_bin(
    N = 134,
    alpha = 0.05,
    se.c = 0.81, se.e = 0.90,
    sp.c = 0.66, sp.e = 0.80,
    psi.d = 0.09,
    psi.nd = 0.14,
    prev = 0.47,
    delta.sp=0.05,
    type.sp = "ni"
  )

  expect_s3_class(pow, "dx_bin")
  expect_true(pow$overall_power >= 0)
  expect_true(pow$overall_power <= 1)

  expect_true(all(c("power_se","power_sp","N_se","N_sp","N_tot") %in% names(pow)))

  expect_equal(pow$hypothesis$specificity, "ni")
})


test_that("power increases with larger N", {

  p_small <- paired_power_opt_bin(
    N = 80,
    alpha = 0.05,
    se.c = 0.81, se.e = 0.90,
    sp.c = 0.66, sp.e = 0.80,
    psi.d = 0.09,
    psi.nd = 0.14,
    prev = 0.47
  )

  p_large <- paired_power_opt_bin(
    N = 250,
    alpha = 0.05,
    se.c = 0.81, se.e = 0.90,
    sp.c = 0.66, sp.e = 0.80,
    psi.d = 0.09,
    psi.nd = 0.14,
    prev = 0.47
  )

  expect_true(p_large$overall_power > p_small$overall_power)
})



test_that("NI specification changes power", {

  p_sup <- paired_power_opt_bin(
    N = 134,
    alpha = 0.05,
    se.c = 0.81, se.e = 0.90,
    sp.c = 0.66, sp.e = 0.80,
    psi.d = 0.09,
    psi.nd = 0.14,
    prev = 0.47
  )

  p_ni <- paired_power_opt_bin(
    N = 134,
    alpha = 0.05,
    se.c = 0.81, se.e = 0.90,
    sp.c = 0.66, sp.e = 0.80,
    psi.d = 0.09,
    psi.nd = 0.14,
    prev = 0.47,
    delta.sp = 0.05,
    type.sp = "ni"
  )

  expect_true(abs(p_sup$overall_power - p_ni$overall_power) > 1e-8)
})

