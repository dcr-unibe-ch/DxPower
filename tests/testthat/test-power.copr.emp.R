test_that("Wald rejection rule equals the Wilson lower limit from prop.test", {
  n <- 37; p0 <- 0.65; alpha <- 0.025
  x <- 0:n
  z <- qnorm(1 - alpha)
  mine <- (x / n - p0) / sqrt(p0 * (1 - p0) / n) > z
  wilson <- vapply(x, function(k) {
    prop.test(k, n, conf.level = 1 - 2 * alpha, correct = FALSE)$conf.int[1] > p0
  }, logical(1))
  expect_identical(mine, wilson)
})

test_that("empirical power approaches the asymptotic power for large n", {
  skip_on_cran()
  set.seed(11)
  B <- 50000
  mc <- function(p) 4 * sqrt(p * (1 - p) / B)
  for (ci in c("wald", "logit")) {
    for (alt in c("one.sided", "two.sided")) {
      a <- list(se.e = 0.85, se.min = 0.80, sp.e = 0.85, sp.min = 0.80,
                ci_method = ci, alternative = alt)
      N <- 1400  # prev 0.5 -> 700 per subgroup
      an  <- do.call(power_opt_bin, c(list(N = N, prev = 0.5), a))
      emp <- do.call(power.copr.emp, c(list(n_dis = 700, n_nondis = 700, n_sim = B), a))
      expect_lt(abs(emp[["se"]] - an$power_se), mc(an$power_se) + 0.01)
      expect_lt(abs(emp[["sp"]] - an$power_sp), mc(an$power_sp) + 0.01)
      expect_lt(abs(emp[["overall"]] - an$power), mc(an$power) + 0.01)
    }
  }
})

test_that("unpaired empirical power approaches the asymptotic power", {
  skip_on_cran()
  set.seed(12)
  B <- 50000
  mc <- function(p) 4 * sqrt(p * (1 - p) / B)
  a <- list(se.c = 0.75, se.e = 0.85, sp.c = 0.65, sp.e = 0.75, alpha = 0.025)
  for (v in c("unpooled", "pooled")) {
    an  <- do.call(power.copr.comp.test,
                   c(list(N = 1000, prev = 0.6, variance = v), a))
    emp <- do.call(power.copr.comp.emp,
                   c(list(n_dis = 600, n_nondis = 400, variance = v, n_sim = B), a))
    expect_lt(abs(emp[["overall"]] - an$power), mc(an$power) + 0.01)
  }
})

test_that("single-test and unpaired empirical power are not the same", {
  skip_on_cran()
  set.seed(13)
  single <- power.copr.emp(100, 100, se.e = 0.85, se.min = 0.75,
                                sp.e = 0.85, sp.min = 0.75, alpha = 0.025,
                                n_sim = 50000)
  unp <- power.copr.comp.emp(100, 100, se.c = 0.75, se.e = 0.85,
                                      sp.c = 0.75, sp.e = 0.85, alpha = 0.025,
                                      n_sim = 50000)
  # the estimated comparator adds variance, so power is clearly lower
  expect_gt(single[["overall"]] - unp[["overall"]], 0.1)
})

test_that("output is a named vector in [0, 1] and reproducible", {
  set.seed(1)
  a <- power.copr.emp(50, 50, 0.85, 0.65, 0.85, 0.55, n_sim = 1000)
  set.seed(1)
  b <- power.copr.emp(50, 50, 0.85, 0.65, 0.85, 0.55, n_sim = 1000)
  expect_equal(a, b)
  expect_named(a, c("overall", "se", "sp"))
  expect_true(all(a >= 0 & a <= 1))
  expect_error(power.copr.emp(50, 50, 0.85, 0.65, 0.85, 0.55, alpha = 2))
})
