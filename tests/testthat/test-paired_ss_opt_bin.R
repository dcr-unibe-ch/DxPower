test_that("paired_ss_opt_bin returns sample size for valid structure for superiority", {

  res <- paired_ss_opt_bin(
    alpha = 0.05,
    power = 0.80,
    se.c = 0.80, se.e = 0.90,
    sp.c = 0.70, sp.e = 0.80,
    psi.d = 0.10,
    psi.nd = 0.15,
    prev = 0.5
  )


  expect_s3_class(res, "dx_bin")
  expect_true(res$N_tot > 0)
  expect_true(res$N_se > 0)
  expect_true(res$N_sp > 0)
  # overall power close to requested
  expect_equal(res$target_power, 0.80, tolerance = 0.02)
  # power is product of components
  expect_equal(res$target_power,
               res$power_se * res$power_sp,
               tolerance = 1e-6)

  expect_true(all(c("power_se","power_sp","N_se","N_sp","N_tot") %in% names(res)))
})


test_that("paired_ss_opt_bin works for non-inferiority", {

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

  expect_equal(res$hypothesis$specificity, "ni")
})




test_that("NI changes required sample size", {

  ss_sup <- paired_ss_opt_bin(
    alpha = 0.05,
    power = 0.80,
    se.c = 0.80, se.e = 0.90,
    sp.c = 0.70, sp.e = 0.75,
    psi.d = 0.10,
    psi.nd = 0.20,
    prev = 0.5
  )

  ss_ni <- paired_ss_opt_bin(
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

  expect_false(identical(ss_sup$N_tot, ss_ni$N_tot))
})



test_that("stricter NI margin increases sample size", {

  ss_loose <- paired_ss_opt_bin(
    alpha = 0.05,
    power = 0.80,
    se.c = 0.80, se.e = 0.90,
    sp.c = 0.66, sp.e = 0.80,
    psi.d = 0.10,
    psi.nd = 0.14,
    prev = 0.47,
    delta.sp = 0.10,
    type.sp = "ni"
  )

  ss_strict <- paired_ss_opt_bin(
    alpha = 0.05,
    power = 0.80,
    se.c = 0.80, se.e = 0.90,
    sp.c = 0.66, sp.e = 0.80,
    psi.d = 0.10,
    psi.nd = 0.14,
    prev = 0.47,
    delta.sp = 0.02,
    type.sp = "ni"
  )

  expect_true(ss_strict$N_tot > ss_loose$N_tot)
})
