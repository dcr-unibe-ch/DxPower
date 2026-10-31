a <- list(se.e = 0.90, se.min = 0.81, sp.e = 0.80, sp.min = 0.66)

# brute-force exact power by enumerating all counts
enum_power <- function(n, e, p0, z, ci) {
  x <- 0:n
  st <- if (ci == "wald") (x / n - p0) / sqrt(p0 * (1 - p0) / n)
        else (qlogis(x / n) - qlogis(p0)) * sqrt(n * p0 * (1 - p0))
  sum(dbinom(x, n, e)[st > z])
}

test_that("exact endpoint power equals enumeration over all counts", {
  for (ci in c("wald", "logit")) {
    for (n in c(1, 5, 15, 40, 129, 300)) {
      z <- qnorm(0.975)
      expect_equal(exact_endpoint_power(n, 0.9, 0.81, z, ci),
                   enum_power(n, 0.9, 0.81, z, ci), tolerance = 1e-10)
    }
  }
  expect_equal(exact_endpoint_power(0, 0.9, 0.81, 1.96), 0)
})

test_that("power_opt_bin(exact = TRUE) is the product of exact endpoint powers", {
  for (ci in c("wald", "logit")) {
    r <- do.call(power_opt_bin, c(list(N = 200, prev = 0.4, ci_method = ci, exact = TRUE), a))
    expect_true(r$exact)
    expect_equal(r$n_dis + r$n_nondis, 200)
    expect_equal(r$power, r$power_se * r$power_sp)
    expect_equal(r$power_se, enum_power(r$n_dis, 0.9, 0.81, qnorm(0.975), ci))
  }
  expect_error(do.call(power_opt_bin, c(list(N = 200.5, prev = 0.4, exact = TRUE), a)),
               "integer")
})

test_that("exact power agrees with empirical power", {
  skip_on_cran()
  set.seed(5)
  B <- 50000
  ex  <- do.call(power_opt_bin, c(list(N = 220, prev = 0.45, exact = TRUE), a))
  emp <- do.call(power.copr.emp,
                 c(list(n_dis = ex$n_dis, n_nondis = ex$n_nondis, n_sim = B), a))
  expect_lt(abs(emp[["overall"]] - ex$power),
            4 * sqrt(ex$power * (1 - ex$power) / B))
})

test_that("both rules match a brute-force scan", {
  for (p in c(0.1, 0.3, 0.5, 0.8)) {
    asy <- do.call(ss_opt_bin, c(list(power = 0.8, prev = p), a))
    Ns <- 2:(3 * asy$n)
    pw <- vapply(Ns, function(N)
      do.call(power_opt_bin, c(list(N = N, prev = p, exact = TRUE), a))$power,
      numeric(1))
    first <- do.call(ss_opt_bin, c(list(power = 0.8, prev = p, exact = TRUE), a))
    stay  <- do.call(ss_opt_bin, c(list(power = 0.8, prev = p, exact = TRUE,
                                        exact_rule = "stays_above"), a))
    expect_equal(first$n, Ns[which(pw >= 0.8)[1]])
    expect_equal(stay$n,  Ns[max(which(pw < 0.8)) + 1])
    expect_gte(stay$n, first$n)
    expect_gte(first$power, 0.8)
    expect_gte(stay$power, 0.8)
    expect_equal(stay$exact_rule, "stays_above")
    expect_equal(stay$target_power, 0.8)
  }
})

test_that("exact_n extends the window when the start value is far off", {
  pw <- function(N) {
    ex <- exact_power_bin(N, 0.4, 0.9, 0.81, 0.8, 0.66, qnorm(0.975), "wald")
    ex$power_se * ex$power_sp
  }
  ref_first <- {Ns <- 2:3000; p <- pw(Ns); Ns[which(p >= 0.8)[1]]}
  ref_stay  <- {Ns <- 2:3000; p <- pw(Ns); Ns[max(which(p < 0.8)) + 1]}
  for (start in c(30, 100, 1500)) {
    expect_equal(exact_n(pw, 0.8, start, "first"), ref_first)
    expect_equal(exact_n(pw, 0.8, start, "stays_above"), ref_stay)
  }
})

test_that("power.copr.test passes exact through", {
  s <- do.call(power.copr.test, c(list(power = 0.8, prev = 0.5, exact = TRUE,
                                       exact_rule = "stays_above"), a))
  expect_equal(s$exact_rule, "stays_above")
  p <- do.call(power.copr.test, c(list(N = s$n, prev = 0.5, exact = TRUE), a))
  expect_equal(p$power, s$power)
})
