test_that("power_opt_bin returns correct structure", {

  pow <- power_opt_bin(
    N=133,
    se.e = 0.85,
    se.min = 0.70,
    sp.e = 0.85,
    sp.min = 0.65,
    alpha = 0.025,
    prev = 0.5
  )

  expect_s3_class(pow, "dx_bin")
  expect_true(all(c("overall_power","power_se","power_sp") %in% names(pow)))
})

