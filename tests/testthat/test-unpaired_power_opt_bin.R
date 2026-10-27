args <- list(N = 910, prev = 0.6, alpha = 0.05,
             se.c = 0.75, se.e = 0.85, sp.c = 0.65, sp.e = 0.75)

# empirical power of the Wald test for a difference of two independent
# proportions, simulated from binomial data
emp_power <- function(n, c, e, margin = 0, B = 20000, pooled = FALSE) {
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

test_that("rho = 0 reproduces the product of endpoint powers", {
  pow <- do.call(unpaired_power_opt_bin, modifyList(args, list(N = 300)))
  expect_equal(pow$overall_power, pow$power_se * pow$power_sp, tolerance = 1e-8)
})

test_that("positive correlation raises and negative lowers joint power", {
  a <- modifyList(args, list(N = 300))
  p0 <- do.call(unpaired_power_opt_bin, a)$overall_power
  pp <- do.call(unpaired_power_opt_bin, c(a, rho = 0.5))$overall_power
  pn <- do.call(unpaired_power_opt_bin, c(a, rho = -0.5))$overall_power
  expect_gt(pp, p0)
  expect_lt(pn, p0)
})

test_that("invalid rho is rejected", {
  expect_error(do.call(unpaired_power_opt_bin, c(args, rho = 1)))
})

test_that("power is consistent with unpaired_ss_opt_bin", {
  res <- unpaired_ss_opt_bin(power = 0.8, prev = 0.6, alpha = 0.05,
                             se.c = 0.75, se.e = 0.85,
                             sp.c = 0.65, sp.e = 0.75)
  pow <- unpaired_power_opt_bin(N = res$N_tot, prev = 0.6, alpha = 0.05,
                                se.c = 0.75, se.e = 0.85,
                                sp.c = 0.65, sp.e = 0.75)
  expect_gte(pow$overall_power, res$target_power)
  expect_lt(pow$overall_power - res$target_power, 0.01)
})

test_that("analytic power matches empirical power (superiority and NI)", {
  skip_on_cran()
  set.seed(2024)
  B <- 20000
  a <- modifyList(args, list(N = 300, delta.sp = 0.05, type.sp = "ni"))
  nd <- a$N * a$prev; nh <- a$N - nd
  for (v in c("pooled", "unpooled")) {
  pow <- do.call(unpaired_power_opt_bin, c(a, variance = v))
  rej_se <- emp_power(nd, a$se.c, a$se.e, B = B, pooled = v == "pooled")
  rej_sp <- emp_power(nh, a$sp.c, a$sp.e, margin = 0.05, B = B,
                      pooled = v == "pooled")
  mc <- function(p) 4 * sqrt(p * (1 - p) / B)
  expect_lt(abs(mean(rej_se) - pow$power_se), mc(pow$power_se))
  expect_lt(abs(mean(rej_sp) - pow$power_sp), mc(pow$power_sp))
  # independent subgroups: joint rejection rate matches rho = 0 power
  expect_lt(abs(mean(rej_se & rej_sp) - pow$overall_power),
            mc(pow$overall_power))
  }
})
