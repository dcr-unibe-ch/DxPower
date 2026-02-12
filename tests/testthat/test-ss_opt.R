test_that("ss_opt returns correct structure", {

  pow <- power_diag(
    N=133,
    se1 = 0.85,
    se0 = 0.70,
    sp1 = 0.85,
    sp0 = 0.65,
    alpha = 0.025,
    prev = 0.5
  )

  expect_type(pow, "list")
  expect_true(all(c("overall","pwr_se","pwr_sp") %in% names(pow)))
})

