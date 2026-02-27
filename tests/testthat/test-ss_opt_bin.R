test_that("ss_opt_bin returns correct structure", {

  res <- ss_opt_bin(
    se.e = 0.85,
    se.min = 0.70,
    sp.e = 0.85,
    sp.min = 0.65,
    alpha = 0.025,
    prev = 0.5,
    power=0.8
  )

  expect_s3_class(res, "dx_bin")
  expect_true(all(c("power_se","power_sp","N_se","N_sp","N_tot") %in% names(res)))
})
