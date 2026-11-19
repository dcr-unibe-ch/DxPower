test_that("Wald rejection rule equals the Wilson lower limit from prop.test", {
  n <- 37; p0 <- 0.65; alpha <- 0.025
  x <- 0:n
  mine <- reject_single(x, n, p0, qnorm(1 - alpha), "wald")
  wilson <- vapply(x, function(k) {
    prop.test(k, n, conf.level = 1 - 2 * alpha, correct = FALSE)$conf.int[1] > p0
  }, logical(1))
  expect_identical(mine, wilson)
})

test_that("reject_single rejects iff x reaches the critical count (exact rule)", {
  z <- qnorm(0.975)
  for (ci in c("wald", "logit")) {
    for (n in c(20, 61, 150)) {
      x <- 0:n
      rej <- reject_single(x, n, 0.7, z, ci)
      # rejection region is an upper tail ...
      expect_false(is.unsorted(rej))
      # ... and its binomial probability is the exact power
      expect_equal(sum(dbinom(x, n, 0.85)[rej]),
                   exact_endpoint_power(n, 0.85, 0.7, z, ci))
    }
  }
})

test_that("reject_single handles the boundary counts", {
  z <- qnorm(0.975)
  expect_false(reject_single(0, 50, 0.7, z, "logit"))   # logit(0) = -Inf
  expect_true(reject_single(50, 50, 0.7, z, "logit"))   # logit(1) = Inf
  expect_false(reject_single(0, 50, 0.7, z, "wald"))
  expect_true(reject_single(50, 50, 0.7, z, "wald"))
})

test_that("reject_unpaired matches hand-computed statistics", {
  # x = 60, y = 70 of n = 100: unpooled z = 0.1 / sqrt(0.45 / 100) = 1.4907,
  # pooled z = 0.1 / sqrt(2 * 0.65 * 0.35 / 100) = 1.4825
  expect_true(reject_unpaired(60, 70, 100, 0, 1.49, "unpooled"))
  expect_false(reject_unpaired(60, 70, 100, 0, 1.49, "pooled"))
  expect_true(reject_unpaired(60, 70, 100, 0, 1.48, "pooled"))
  # non-inferiority: x = 70, y = 65, margin 0.1: z = 0.05 / sqrt(0.4375 / 100) = 0.756
  expect_true(reject_unpaired(70, 65, 100, 0.1, 0.75, "unpooled"))
  expect_false(reject_unpaired(70, 65, 100, 0.1, 0.76, "unpooled"))
  expect_false(reject_unpaired(70, 65, 100, 0, 0.75, "unpooled"))  # no margin
})

test_that("reject_unpaired treats a zero standard error as not rejected", {
  for (v in c("unpooled", "pooled")) {
    expect_false(reject_unpaired(0, 0, 50, 0, 1.96, v))
    expect_false(reject_unpaired(50, 50, 50, 0, 1.96, v))
  }
  # vectorised
  expect_identical(reject_unpaired(c(0, 10), c(0, 40), 50, 0, 1.96, "unpooled"),
                   c(FALSE, TRUE))
})

test_that("reject_unpaired reproduces the exact unpaired power", {
  n <- 30; alpha <- 0.025; z <- qnorm(1 - alpha)
  g <- expand.grid(x = 0:n, y = 0:n)
  for (v in c("unpooled", "pooled")) {
    for (type in c("sup", "ni")) {
      margin <- if (type == "ni") 0.1 else 0
      rej <- reject_unpaired(g$x, g$y, n, margin, z, v)
      p <- sum(dbinom(g$x, n, 0.7) * dbinom(g$y, n, 0.8) * rej)
      expect_equal(p, unpaired_exact_endpoint_power(n, alpha, 0.7, 0.8, 0.1,
                                                    type, v))
    }
  }
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
