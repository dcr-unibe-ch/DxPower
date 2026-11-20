args <- list(prev = 0.5, se.e = 0.85, se.min = 0.65,
             sp.e = 0.85, sp.min = 0.55)

test_that("power given: matches ss_opt_bin()", {
  a <- do.call(power.copr.test, c(list(power = 0.8), args))
  b <- do.call(ss_opt_bin, c(list(power = 0.8), args))
  expect_equal(a, b)
  expect_equal(a$target_power, 0.8)
})

test_that("N given: matches power_opt_bin()", {
  a <- do.call(power.copr.test, c(list(N = 122, ci_method = "logit"), args))
  b <- do.call(power_opt_bin, c(list(N = 122, ci_method = "logit"), args))
  expect_equal(a, b)
  expect_null(a$target_power)
})

test_that("round trip: power at the computed N is about the target", {
  ss <- do.call(power.copr.test, c(list(power = 0.8), args))
  pw <- do.call(power.copr.test, c(list(N = ceiling(ss$n)), args))
  expect_gte(pw$power, 0.8 - 1e-3)
  expect_lt(pw$power, 0.85)
})

test_that("exactly one of N and power must be NULL", {
  expect_error(do.call(power.copr.test, args), "exactly one")
  expect_error(do.call(power.copr.test, c(list(N = 100, power = 0.8), args)),
               "exactly one")
})
