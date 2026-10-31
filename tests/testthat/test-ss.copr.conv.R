a <- list(prev = 0.2, se.e = 0.85, se.min = 0.65, sp.e = 0.85, sp.min = 0.55)

test_that("N is the maximum of the two endpoint requirements", {
  res <- do.call(ss.copr.conv, c(list(power = 0.8), a))
  expect_s3_class(res, "power.htest")
  expect_equal(res$n, max(res$N_se_required, res$N_sp_required))
  expect_match(res$method, "conventional")
  expect_equal(res$power_endpoint, 0.9)
})

test_that("each endpoint reaches the endpoint power at its required N", {
  for (p in c(0.2, 0.5, 0.8)) {
    for (ci in c("wald", "logit")) {
      args <- modifyList(a, list(prev = p))
      res <- do.call(ss.copr.conv, c(list(power = 0.8, ci_method = ci), args))
      pw_se <- do.call(power_opt_bin, c(list(N = res$N_se_required, ci_method = ci), args))
      pw_sp <- do.call(power_opt_bin, c(list(N = res$N_sp_required, ci_method = ci), args))
      expect_gte(pw_se$power_se, 0.9 - 1e-12)
      expect_gte(pw_sp$power_sp, 0.9 - 1e-12)
    }
  }
})

test_that("conventional design is never smaller than, and overpowers, the optimal one", {
  for (p in c(0.1, 0.3, 0.5, 0.7, 0.9)) {
    args <- modifyList(a, list(prev = p))
    conv <- do.call(ss.copr.conv, c(list(power = 0.8), args))
    opt  <- do.call(ss_opt_bin,  c(list(power = 0.8), args))
    expect_gte(conv$n, opt$n)
    expect_gte(conv$power, 0.8)
  }
})

test_that("invalid power_endpoint is rejected", {
  expect_error(do.call(ss.copr.conv, c(list(power = 0.8, power_endpoint = 0.7), a)),
               "power_endpoint")
})
