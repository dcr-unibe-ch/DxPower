un <- list(prev = 0.6, alpha = 0.025, se.c = 0.75, se.e = 0.85,
           sp.c = 0.65, sp.e = 0.75)
pa <- list(prev = 0.47, alpha = 0.025, se.c = 0.81, se.e = 0.90,
           sp.c = 0.66, sp.e = 0.80, psi.d = 0.09, psi.nd = 0.14,
           design = "paired")

pwr_at <- function(N, ...) {
  do.call(power.copr.comp.test, c(un, list(N = N), list(...)))$power
}

# -- class and method ---------------------------------------------------------

test_that("results are power.htest objects whose method states the design", {
  u <- do.call(power.copr.comp.test, c(un, list(N = 300)))
  p <- do.call(power.copr.comp.test, c(pa, list(N = 134)))
  expect_s3_class(u, "power.htest")
  expect_s3_class(p, "power.htest")
  expect_match(u$method, "Comparative unpaired design")
  expect_match(p$method, "Comparative paired design")
  expect_match(u$method, "normal approximation")
  expect_match(p$method, "normal approximation")
  expect_match(u$method, "Wald|Logit")
  expect_match(u$method, "power calculation")
  expect_match(do.call(power.copr.comp.test, c(un, list(power = 0.8)))$method,
               "sample size calculation")
  expect_output(print(u), "Comparative unpaired design")
})

test_that("all single-test functions return power.htest with a descriptive method", {
  a <- list(prev = 0.5, se.e = 0.85, se.min = 0.65, sp.e = 0.85, sp.min = 0.55)
  res <- list(
    do.call(power_opt_bin, c(list(N = 122), a)),
    do.call(power_opt_bin, c(list(N = 122, exact = TRUE), a)),
    do.call(ss_opt_bin, c(list(power = 0.8), a)),
    do.call(ss_opt_bin, c(list(power = 0.8, exact = TRUE), a)),
    do.call(ss.copr.conv, c(list(power = 0.8), a)),
    do.call(power.copr.test, c(list(power = 0.8), a))
  )
  for (r in res) {
    expect_s3_class(r, "power.htest")
    expect_match(r$method, "Single-test design")
    expect_match(r$method, "Wald|Logit")
    expect_match(r$method, "normal approximation|Exact")
    expect_true(all(c("n", "power", "sig.level", "power_se", "power_sp") %in% names(r)))
  }
  expect_match(res[[1]]$method, "normal approximation")
  expect_match(res[[2]]$method, "Exact")
  expect_match(res[[3]]$method, "sample size")
  expect_match(res[[5]]$method, "conventional")
})

# -- argument handling --------------------------------------------------------

test_that("exactly one of N and power must be NULL", {
  expect_error(do.call(power.copr.comp.test, un), "exactly one")
  expect_error(do.call(power.copr.comp.test, c(un, list(N = 100, power = 0.8))),
               "exactly one")
})

test_that("paired design defaults to the lower Stark bound with a message", {
  p <- pa[setdiff(names(pa), c("psi.d", "psi.nd"))]
  expect_message(r <- do.call(power.copr.comp.test, c(p, list(N = 134))),
                 "lower bound of Stark")
  # |se.c - se.e| = 0.09 and |sp.c - sp.e| = 0.14 are the psi of `pa`
  expect_equal(r$psi.d, 0.09)
  expect_equal(r$psi.nd, 0.14)
  expect_equal(r$power, do.call(power.copr.comp.test, c(pa, list(N = 134)))$power)
  # only the missing endpoint is defaulted
  expect_message(r2 <- do.call(power.copr.comp.test,
                               c(p, list(N = 134, psi.d = 0.12))),
                 "psi.nd and rho.nd not given")
  expect_equal(r2$psi.d, 0.12)
})

test_that("psi is checked against its bounds and ignored (with warning) if unpaired", {
  expect_error(do.call(power.copr.comp.test,
                       c(modifyList(pa, list(psi.d = 0.5)), list(N = 134))),
               "psi.d outside")
  expect_warning(do.call(power.copr.comp.test, c(un, list(N = 300, psi.d = 0.1))),
                 "only used for design")
  expect_warning(do.call(power.copr.comp.test,
                         c(pa, list(N = 134, variance = "pooled"))),
                 "only used for design")
})

# -- paired: correlation instead of discordance -------------------------------

test_that("rho_to_psi and psi_to_rho are inverses within the Stark range", {
  for (rho in c(0, 0.3, 0.6)) {
    psi <- rho_to_psi(0.81, 0.90, rho)
    expect_equal(psi_to_rho(0.81, 0.90, psi), rho)
  }
  # independence: psi = p1 + p2 - 2 p1 p2
  expect_equal(rho_to_psi(0.81, 0.90, 0), 0.81 + 0.90 - 2 * 0.81 * 0.90)
  b <- paired_bounds(0.81, 0.90, 0.66, 0.80)
  expect_equal(names(b), c("endpoint", "psi.min", "psi.max", "rho.min", "rho.max"))
  # Stark et al. (2022, Eq. 10): |c - e| <= psi <= c + e - 2ce, i.e. rho in [0, rho.max]
  expect_equal(b$psi.min, c(0.09, 0.14))
  expect_equal(b$psi.max, c(0.81 + 0.90 - 2 * 0.81 * 0.90,
                            0.66 + 0.80 - 2 * 0.66 * 0.80))
  expect_equal(b$rho.min, c(0, 0))
  expect_equal(rho_to_psi(0.81, 0.90, b$rho.max[1]), b$psi.min[1])
  expect_equal(rho_to_psi(0.81, 0.90, b$rho.min[1]), b$psi.max[1])
  expect_equal(psi_to_rho(0.81, 0.90, b$psi.min[1]), b$rho.max[1])
  expect_equal(psi_to_rho(0.81, 0.90, b$psi.max[1]), 0)
  # outside the range (including negative correlation) is an error
  expect_error(rho_to_psi(0.81, 0.90, b$rho.max[1] + 0.01), "range of Stark")
  expect_error(rho_to_psi(0.81, 0.90, -0.05), "range of Stark")
  expect_error(psi_to_rho(0.81, 0.90, 0.05), "range of Stark")
  expect_error(psi_to_rho(0.81, 0.90, b$psi.max[1] + 0.01), "range of Stark")
})

test_that("paired: psi above c + e - 2ce (negative correlation) is an error", {
  base <- pa[setdiff(names(pa), c("psi.d", "psi.nd"))]
  upper <- paired_bounds(0.81, 0.90, 0.66, 0.80)$psi.max[1]
  expect_error(do.call(power.copr.comp.test,
                       c(base, list(psi.d = upper + 0.01, psi.nd = 0.14, N = 134))),
               "psi.d outside")
  expect_error(do.call(power.copr.comp.test,
                       c(base, list(rho.d = -0.05, psi.nd = 0.14, N = 134))),
               "range of Stark")
  expect_no_warning(do.call(power.copr.comp.test,
                            c(base, list(psi.d = upper, psi.nd = 0.14, N = 134))))
})

test_that("paired: rho gives the same result as the equivalent psi", {
  rho.d <- 0.5; rho.nd <- 0.3
  psi.d <- rho_to_psi(pa$se.c, pa$se.e, rho.d)
  psi.nd <- rho_to_psi(pa$sp.c, pa$sp.e, rho.nd)
  base <- pa[setdiff(names(pa), c("psi.d", "psi.nd"))]
  for (extra in list(list(), list(delta.sp = 0.05, type.sp = "ni"))) {
    r <- do.call(power.copr.comp.test,
                 c(base, extra, list(rho.d = rho.d, rho.nd = rho.nd, N = 134)))
    p <- do.call(power.copr.comp.test,
                 c(base, extra, list(psi.d = psi.d, psi.nd = psi.nd, N = 134)))
    expect_equal(r$power, p$power)
    expect_equal(r$psi.d, psi.d)
    expect_equal(r$psi.nd, psi.nd)
  }
  # sample size too, and psi / rho may be mixed across endpoints
  expect_equal(
    do.call(power.copr.comp.test,
            c(base, list(rho.d = rho.d, psi.nd = psi.nd, power = 0.8)))$n,
    do.call(power.copr.comp.test,
            c(base, list(psi.d = psi.d, psi.nd = psi.nd, power = 0.8)))$n)
  # larger positive correlation (smaller discordance) needs fewer subjects
  n_rho <- function(r) do.call(power.copr.comp.test,
                               c(base, list(rho.d = r, rho.nd = r, power = 0.8)))$n
  expect_lt(n_rho(0.6), n_rho(0.3))
})

test_that("paired: exactly one of psi and rho per endpoint", {
  base <- pa[setdiff(names(pa), c("psi.d", "psi.nd"))]
  expect_error(do.call(power.copr.comp.test,
                       c(base, list(psi.d = 0.09, rho.d = 0.5, rho.nd = 0.3, N = 134))),
               "only one of psi.d and rho.d")
  # the other endpoint falls back to the lower Stark bound
  expect_message(do.call(power.copr.comp.test, c(base, list(rho.d = 0.5, N = 134))),
                 "psi.nd and rho.nd not given")
  expect_error(do.call(power.copr.comp.test,
                       c(base, list(rho.d = 0.99, rho.nd = 0.3, N = 134))),
               "rho.d = 0.99 is outside the range of Stark")
  expect_warning(do.call(power.copr.comp.test, c(un, list(N = 300, rho.d = 0.1))),
                 "only used for design")
})

# -- unpaired -----------------------------------------------------------------

test_that("unpaired: overall power is the product of endpoint powers", {
  pow <- do.call(power.copr.comp.test, c(un, list(N = 300)))
  expect_equal(pow$power, pow$power_se * pow$power_sp, tolerance = 1e-8)
})

test_that("unpaired: n is the smallest N reaching the target power", {
  for (extra in list(list(), list(delta.sp = 0.05, type.sp = "ni"))) {
    res <- do.call(power.copr.comp.test, c(un, list(power = 0.8), extra))
    expect_gte(do.call(pwr_at, c(list(res$n), extra)), 0.8)
    expect_lt(do.call(pwr_at, c(list(res$n - 1), extra)), 0.8)
    expect_equal(res$target_power, 0.8)
    expect_gte(res$power, res$target_power)
    expect_lt(res$power - res$target_power, 0.01)
    expect_equal(res$n_dis + res$n_nondis, res$n)
  }
})

test_that("unpaired: power is non-decreasing in N", {
  p <- vapply(seq(50, 1500, by = 25), pwr_at, numeric(1))
  expect_true(all(diff(p) >= -1e-10))
})

test_that("find_n errors when the target is unreachable", {
  expect_error(find_n(function(n) 0.5, 0.8), "Upper bound")
})

test_that("unpaired: analytic power matches empirical power (sup and NI)", {
  skip_on_cran()
  set.seed(2024)
  B <- 20000
  # empirical power of the Wald test for a difference of two independent
  # proportions, simulated from binomial data
  emp_power <- function(n, c, e, margin = 0, pooled = FALSE) {
    x <- rbinom(B, n, c) / n
    y <- rbinom(B, n, e) / n
    se <- if (pooled) {
      pb <- (x + y) / 2
      sqrt(2 * pb * (1 - pb) / n)
    } else {
      sqrt((x * (1 - x) + y * (1 - y)) / n)
    }
    (y - x + margin) / se > qnorm(1 - 0.025)
  }
  mc <- function(p) 4 * sqrt(p * (1 - p) / B)
  a <- c(un, list(N = 300, delta.sp = 0.05, type.sp = "ni"))
  nd <- a$N * a$prev; nh <- a$N - nd
  for (v in c("pooled", "unpooled")) {
    pow <- do.call(power.copr.comp.test, c(a, variance = v))
    rej_se <- emp_power(nd, a$se.c, a$se.e, pooled = v == "pooled")
    rej_sp <- emp_power(nh, a$sp.c, a$sp.e, margin = 0.05, pooled = v == "pooled")
    expect_lt(abs(mean(rej_se) - pow$power_se), mc(pow$power_se))
    expect_lt(abs(mean(rej_sp) - pow$power_sp), mc(pow$power_sp))
    # independent subgroups: joint rejection rate matches the product of powers
    expect_lt(abs(mean(rej_se & rej_sp) - pow$power), mc(pow$power))
  }
})

# -- paired -------------------------------------------------------------------

test_that("paired: power is in [0, 1] and has the expected structure", {
  for (extra in list(list(), list(delta.sp = 0.05, type.sp = "ni"))) {
    pow <- do.call(power.copr.comp.test, c(pa, list(N = 134), extra))
    expect_gte(pow$power, 0)
    expect_lte(pow$power, 1)
    expect_equal(pow$power, pow$power_se * pow$power_sp)
    expect_true(all(c("power_se", "power_sp", "n_dis", "n_nondis", "n") %in% names(pow)))
  }
  expect_equal(do.call(power.copr.comp.test,
                       c(pa, list(N = 134, delta.sp = 0.05, type.sp = "ni")))$type.sp,
               "ni")
})

test_that("paired: power increases with N and NI changes power", {
  p <- function(N, ...) do.call(power.copr.comp.test, c(pa, list(N = N), list(...)))$power
  expect_gt(p(250), p(80))
  expect_gt(abs(p(134) - p(134, delta.sp = 0.05, type.sp = "ni")), 1e-8)
})

test_that("paired: sample size reaches the target power", {
  for (extra in list(list(), list(delta.sp = 0.05, type.sp = "ni"))) {
    res <- do.call(power.copr.comp.test,
                   c(modifyList(pa, list(prev = 0.5, psi.nd = 0.2, sp.e = 0.75,
                                         se.c = 0.8, psi.d = 0.1)),
                     list(power = 0.8), extra))
    expect_gt(res$n, 0)
    expect_equal(res$target_power, 0.8)
    expect_equal(res$power, 0.8, tolerance = 0.02)
    # round trip: power at the returned N is the target
    back <- do.call(power.copr.comp.test,
                    c(modifyList(pa, list(prev = 0.5, psi.nd = 0.2, sp.e = 0.75,
                                          se.c = 0.8, psi.d = 0.1)),
                      list(N = res$n), extra))
    expect_equal(back$power, res$power)
  }
  res <- do.call(power.copr.comp.test, c(pa, list(power = 0.8)))
  back <- do.call(power.copr.comp.test, c(pa, list(N = res$n)))
  expect_equal(res$n, round(res$n))
  expect_gte(back$power, 0.80)
  expect_lt(back$power, 0.81)
  # n is the smallest integer reaching the target
  expect_lt(do.call(power.copr.comp.test, c(pa, list(N = res$n - 1)))$power, 0.80)
})

test_that("paired: NI and its margin change the required sample size", {
  ss <- function(...) do.call(power.copr.comp.test,
                              c(modifyList(pa, list(psi.nd = 0.14)),
                                list(power = 0.8), list(...)))$n
  expect_false(identical(ss(), ss(delta.sp = 0.05, type.sp = "ni")))
  expect_gt(ss(delta.sp = 0.02, type.sp = "ni"), ss(delta.sp = 0.10, type.sp = "ni"))
})

# -- alternative --------------------------------------------------------------

test_that("two-sided is one-sided at alpha / 2 (both designs)", {
  for (d in list(un, pa)) {
    two <- do.call(power.copr.comp.test,
                   c(modifyList(d, list(alpha = 0.05)),
                     list(N = 200, alternative = "two.sided")))
    one <- do.call(power.copr.comp.test, c(d, list(N = 200)))
    expect_equal(two$power, one$power)
    expect_equal(two$alternative, "two.sided")
    expect_equal(two$sig.level, 0.05)
    expect_equal(one$alternative, "one.sided")
  }
  n2 <- do.call(power.copr.comp.test, c(un, list(power = 0.8, alternative = "two.sided")))$n
  n1 <- do.call(power.copr.comp.test, c(un, list(power = 0.8)))$n
  expect_gt(n2, n1)
})

# -- exact unpaired -----------------------------------------------------------

test_that("exact unpaired: power matches the empirical power", {
  skip_on_cran()
  set.seed(11)
  B <- 1e5
  for (v in c("unpooled", "pooled")) {
    ex <- do.call(power.copr.comp.test, c(un, list(N = 200, exact = TRUE, variance = v)))
    em <- power.copr.comp.emp(ex$n_dis, ex$n_nondis, un$se.c, un$se.e,
                              un$sp.c, un$sp.e, alpha = un$alpha,
                              variance = v, n_sim = B)
    mc <- function(p) 4 * sqrt(p * (1 - p) / B)
    expect_lt(abs(em[["overall"]] - ex$power), mc(ex$power))
    expect_lt(abs(em[["se"]] - ex$power_se), mc(ex$power_se))
    expect_lt(abs(em[["sp"]] - ex$power_sp), mc(ex$power_sp))
  }
})

test_that("exact unpaired: structure, subgroups and NI", {
  ex <- do.call(power.copr.comp.test,
                c(un, list(N = 201, exact = TRUE, delta.sp = 0.05, type.sp = "ni")))
  expect_equal(ex$n_dis, round(201 * un$prev))
  expect_equal(ex$n_dis + ex$n_nondis, 201)
  expect_equal(ex$power, ex$power_se * ex$power_sp)
  expect_match(ex$method, "Exact binomial")
  expect_true(ex$exact)
  expect_equal(ex$variance, "unpooled")
})

test_that("exact unpaired: sample size rules", {
  first <- do.call(power.copr.comp.test, c(un, list(power = 0.8, exact = TRUE)))
  stays <- do.call(power.copr.comp.test,
                   c(un, list(power = 0.8, exact = TRUE, exact_rule = "stays_above")))
  pw <- function(N) do.call(power.copr.comp.test, c(un, list(N = N, exact = TRUE)))$power
  expect_gte(first$power, 0.8)
  expect_lt(pw(first$n - 1), 0.8)
  expect_gte(stays$n, first$n)
  expect_true(all(vapply(stays$n:(stays$n + 40), pw, numeric(1)) >= 0.8))
  expect_lt(pw(stays$n - 1), 0.8)
  expect_equal(first$exact_rule, "first")
})

test_that("exact unpaired: argument checks", {
  expect_error(do.call(power.copr.comp.test, c(un, list(N = 100.5, exact = TRUE))),
               "integer")
  expect_error(do.call(power.copr.comp.test,
                       c(un, list(N = 100, exact = TRUE, variance = "null"))),
               "pooled")
  expect_error(do.call(power.copr.comp.test, c(pa, list(N = 100, exact = TRUE))),
               "unpaired")
})

# -- input validation ---------------------------------------------------------

test_that("input validation is consistent across functions", {
  a <- list(prev = 0.5, se.e = 0.85, se.min = 0.65, sp.e = 0.85, sp.min = 0.55)
  bad <- function(f, ..., msg) expect_error(do.call(f, modifyList(a, list(...))), msg)
  for (f in list(power.copr.test, ss.copr.conv)) {
    bad(f, power = 0.8, prev = 1.2, msg = "prev must be")
    bad(f, power = 0.8, alpha = 0, msg = "alpha must be")
  }
  bad(power.copr.test, N = -5, msg = "N must be")
  bad(power.copr.test, N = c(100, 200), msg = "N must be")
  bad(power.copr.test, N = 100, se.e = 1.1, msg = "se.e must be")
  bad(ss.copr.conv, power = 1, msg = "power must be")
  bad(power.copr.test, N = 100, se.e = 0.5, msg = "se.e must be greater")
  c2 <- c(un, list(N = 100))
  for (nm in c("prev", "alpha", "se.c", "sp.e")) {
    expect_error(do.call(power.copr.comp.test, modifyList(c2, setNames(list(1.5), nm))),
                 paste0(nm, " must be"))
  }
  expect_error(do.call(power.copr.comp.test, modifyList(c2, list(N = 0))), "N must be")
  expect_error(do.call(power.copr.comp.test, modifyList(c2, list(N = NULL, power = 1))),
               "power must be")
})
