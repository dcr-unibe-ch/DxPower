a <- list(prev = 0.6, alpha = 0.05, se.c = 0.75, se.e = 0.85,
          sp.c = 0.65, sp.e = 0.75)
pwr_at <- function(N, ...) {
  do.call(unpaired_power_opt_bin, c(a, list(N = N), list(...)))$overall_power
}

test_that("find_n returns the smallest N reaching the target power", {
  for (extra in list(list(), list(rho = 0.4), list(delta.sp = 0.05, type.sp = "ni"))) {
    res <- do.call(unpaired_ss_opt_bin, c(a, list(power = 0.8), extra))
    expect_s3_class(res, "dx_bin")
    expect_gte(do.call(pwr_at, c(list(res$N_tot), extra)), 0.8)
    expect_lt(do.call(pwr_at, c(list(res$N_tot - 1), extra)), 0.8)
    expect_equal(res$target_power, 0.8)
  }
})

test_that("higher correlation needs a smaller sample", {
  n0 <- do.call(unpaired_ss_opt_bin, c(a, power = 0.8, rho = 0))$N_tot
  n5 <- do.call(unpaired_ss_opt_bin, c(a, power = 0.8, rho = 0.5))$N_tot
  expect_lte(n5, n0)
})

test_that("find_n errors when the target is unreachable", {
  expect_error(find_n(function(n) 0.5, 0.8), "Upper bound")
})

test_that("power is non-decreasing in N", {
  p <- vapply(seq(50, 1500, by = 25), pwr_at, numeric(1))
  expect_true(all(diff(p) >= -1e-10))
})
