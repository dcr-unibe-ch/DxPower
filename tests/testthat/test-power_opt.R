test_that("power_opt returns correct structure", {

  res <- power_opt(
    se1 = 0.85,
    se0 = 0.70,
    sp1 = 0.85,
    sp0 = 0.65,
    alpha = 0.025,
    prev = 0.5,
    pwr=0.8
  )

  expect_type(res, "list")
  expect_true(all(c("pwr_se","pwr_sp","N_se","N_sp","N_tot") %in% names(res)))
})
