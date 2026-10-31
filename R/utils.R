# Smallest integer N with pwr_fun(N) >= power, for a power function that is
# non-decreasing in N: bracket by doubling, then bisect.
find_n <- function(pwr_fun, power, n_low = 2, n_high = 1000, n_max = 1e7) {
  if (pwr_fun(n_low) >= power) {
    warning("Target power is already reached at n_low = ", n_low, ".")
    return(n_low)
  }

  while (pwr_fun(n_high) < power) {
    n_low <- n_high
    n_high <- n_high * 2
    if (n_high > n_max) stop("Upper bound exploded; check inputs.")
  }

  lo <- n_low
  hi <- n_high
  while (lo + 1 < hi) {
    mid <- floor((lo + hi) / 2)
    if (pwr_fun(mid) >= power) hi <- mid else lo <- mid
  }
  hi
}

# Effect size and variances (per observation) of one endpoint on the scale used
# by the test: probability scale (Wald) or logit scale. `v_min` is the variance
# under the null value, `v_e` under the alternative.
endpoint_scale <- function(e, min, ci_method = c("wald", "logit")) {
  ci_method <- match.arg(ci_method)
  if (ci_method == "wald") {
    list(d = e - min, v_min = min * (1 - min), v_e = e * (1 - e))
  } else {
    list(d = qlogis(e) - qlogis(min),
         v_min = 1 / (min * (1 - min)), v_e = 1 / (e * (1 - e)))
  }
}

# Exact binomial power of one endpoint test, vectorised over n.
# The Wald (null variance) and logit statistics are increasing in the count x,
# so the test rejects iff x >= cc with cc the smallest integer such that the
# statistic exceeds z; power = P(X >= cc), X ~ Bin(n, e).
exact_endpoint_power <- function(n, e, min, z, ci_method = c("wald", "logit")) {
  ci_method <- match.arg(ci_method)
  ok <- n >= 1
  n1 <- pmax(n, 1)
  v0 <- min * (1 - min)
  # x must exceed `bound`; the tolerance keeps x == bound from rejecting
  bound <- if (ci_method == "wald") {
    n1 * min + z * sqrt(n1 * v0)
  } else {
    n1 * plogis(qlogis(min) + z / sqrt(n1 * v0))
  }
  cc <- floor(bound + 1e-9) + 1
  pw <- pbinom(cc - 1, n1, e, lower.tail = FALSE)
  pw[!ok] <- 0
  pw
}

# Exact endpoint powers for total sample sizes N (integer, vectorised) split
# into round(N * prev) diseased and the rest non-diseased.
exact_power_bin <- function(N, prev, se.e, se.min, sp.e, sp.min, z, ci_method) {
  n_dis    <- round(N * prev)
  n_nondis <- N - n_dis
  list(n_dis = n_dis, n_nondis = n_nondis,
       power_se = exact_endpoint_power(n_dis,    se.e, se.min, z, ci_method),
       power_sp = exact_endpoint_power(n_nondis, sp.e, sp.min, z, ci_method))
}

# Smallest integer N reaching `power` when the power function `pw` (vectorised
# over N) is not monotone in N. The scan starts in a window around the
# asymptotic sample size N_asym and is extended when needed.
#   rule "first":       first N with pw(N) >= power
#   rule "stays_above": smallest N such that pw(N') >= power for all N' >= N
#                       (within the scanned range)
exact_n <- function(pw, power, N_asym, rule = c("first", "stays_above"),
                    n_max = 1e7) {
  rule <- match.arg(rule)
  lo <- max(2, floor(N_asym / 2))
  hi <- ceiling(1.5 * N_asym) + 10
  p <- pw(lo:hi)

  extend_up <- function() {
    if (hi >= n_max) stop("Upper bound exploded; check inputs.")
    new_hi <- min(2 * hi, n_max)
    p <<- c(p, pw((hi + 1):new_hi))
    hi <<- new_hi
  }
  extend_down <- function() {
    new_lo <- max(2, floor(lo / 2))
    p <<- c(pw(new_lo:(lo - 1)), p)
    lo <<- new_lo
  }

  repeat {
    ok <- p >= power
    if (rule == "first") {
      if (!any(ok)) {
        extend_up()
      } else if (ok[1] && lo > 2) {
        extend_down()
      } else {
        return(lo + which(ok)[1] - 1)
      }
    } else {
      fail <- which(!ok)
      if (!ok[length(ok)]) {
        extend_up()
      } else if (length(fail) == 0) {
        if (lo > 2) extend_down() else return(lo)
      } else {
        return(lo + max(fail))
      }
    }
  }
}

# Result object shared by all power / sample size functions: an object of class
# "power.htest" (printed by stats:::print.power.htest). `n` is the total sample
# size and `power` the achieved overall power; `...` are design-specific
# extras (e.g. ci_method, type.se). NULL entries are dropped.
power_htest <- function(method, n, n_dis, n_nondis, prev, sig.level,
                        alternative, power, power_se, power_sp,
                        target_power = NULL, ...) {
  res <- list(n = n, n_dis = n_dis, n_nondis = n_nondis, prev = prev,
              sig.level = sig.level, alternative = alternative,
              power = power, power_se = power_se, power_sp = power_sp,
              target_power = target_power, ...)
  res <- res[!vapply(res, is.null, logical(1))]
  res$note <- "n is the total sample size."
  res$method <- method
  structure(res, class = "power.htest")
}

# "method" string of a power.htest: design, approximation, co-primary IUT.
power_method <- function(design, approx, analysis = c("power", "samplesize")) {
  analysis <- match.arg(analysis)
  l1 <- paste0(design, " (",
               if (analysis == "power") "power calculation"
               else "sample size calculation", ")")
  # print.power.htest indents only the first line (by 5 characters), so the
  # shorter line is padded to centre the two lines on each other
  w <- nchar(l1) - nchar(approx)
  if (w >= 0) {
    paste0(l1, "\n", strrep(" ", 5 + w %/% 2), approx)
  } else {
    paste0(strrep(" ", (-w) %/% 2), l1, "\n", strrep(" ", 5), approx)
  }
}

# approximation label of the single-test design
single_approx <- function(ci_method, exact) {
  label <- if (ci_method == "wald") "Wald" else "Logit"
  if (exact) paste0("Exact binomial (", label, " test)")
  else paste0(label, " normal approximation")
}
