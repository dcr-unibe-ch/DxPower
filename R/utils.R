#' Smallest sample size reaching a target power
#'
#' Finds the smallest integer total sample size \eqn{N}{N} with
#' \code{pwr_fun(N) >= power}. The upper bound is found by doubling, then the
#' solution by bisection, so \code{pwr_fun} must be non-decreasing in \eqn{N}{N}
#' (for sawtooth power functions see \code{exact_n()}).
#'
#' @param pwr_fun Function of a single sample size \code{N} returning the power
#'   (a number in \eqn{[0, 1]}{[0, 1]}).
#' @param power Target power.
#' @param n_low Lower end of the search interval. A warning is given, and
#'   \code{n_low} returned, if the target is already reached there.
#' @param n_high Initial upper end of the search interval, doubled until the
#'   target is reached.
#' @param n_max Largest upper end tried before stopping with an error (target
#'   unreachable).
#'
#' @return The smallest integer \eqn{N}{N} with \code{pwr_fun(N) >= power}.
#' @keywords internal
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

# -----------------------------------------------------------------------------
# Input validation shared by all power / sample size functions
# -----------------------------------------------------------------------------

# `x` must be a single non-missing number with lower < x < upper.
check_num <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= lower || x >= upper) {
    stop(name, " must be a single number ",
         if (is.finite(lower) && is.finite(upper))
           paste0("between ", lower, " and ", upper)
         else "greater than 0", ".", call. = FALSE)
  }
  invisible(x)
}

# `x` must be a single positive whole number (a count such as n_dis or n_sim).
check_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1 || x != round(x)) {
    stop(name, " must be a single positive integer.", call. = FALSE)
  }
  invisible(x)
}

# Warns if the effective one-sided level exceeds the conventional 0.025
# (alpha = 0.05 with one.sided, or alpha > 0.05 with two.sided). Called once
# from each exported function; silent for an invalid alpha (rejected by the
# validation).
warn_alpha <- function(alpha, alternative = c("one.sided", "two.sided")) {
  alternative <- alternative[1]
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) return(invisible(FALSE))
  alpha_eff <- if (identical(alternative, "two.sided")) alpha / 2 else alpha
  if (alpha_eff > 0.025 + 1e-12) {
    warning("alpha = ", alpha, " with alternative = \"", alternative,
            "\" is a one-sided level of ", alpha_eff, ", above the conventional ",
            "0.025 (alpha = 0.025 one-sided, or 0.05 two-sided). Check that ",
            "this is intended.", call. = FALSE)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

# Exactly one of N and power is NULL (it is the one to calculate).
check_one_of <- function(N, power) {
  if (sum(vapply(list(N, power), is.null, logical(1))) != 1) {
    stop("exactly one of 'N' and 'power' must be NULL", call. = FALSE)
  }
}

# Checks the non-NULL one of N / power, prevalence, alpha and the named
# probabilities.
check_inputs <- function(N = NULL, power = NULL, prev, alpha, probs = list()) {
  if (!is.null(N)) check_num(N, "N", lower = 0)
  if (!is.null(power)) check_num(power, "power", 0, 1)
  check_num(prev, "prev", 0, 1)
  check_num(alpha, "alpha", 0, 1)
  for (nm in names(probs)) check_num(probs[[nm]], nm, 0, 1)
  invisible(TRUE)
}

# Superiority against a fixed null value needs e > min (single-test design).
check_single <- function(se.e, se.min, sp.e, sp.min) {
  for (nm in c("se.e", "se.min", "sp.e", "sp.min")) {
    check_num(get(nm), nm, 0, 1)
  }
  if (se.e <= se.min) stop("se.e must be greater than se.min.", call. = FALSE)
  if (sp.e <= sp.min) stop("sp.e must be greater than sp.min.", call. = FALSE)
  invisible(TRUE)
}

# Critical value of the standard normal test; `alpha` is the one-sided level
# (or the two-sided level, which is halved).
z_critical <- function(alpha, alternative) {
  alpha_eff <- if (alternative == "one.sided") alpha else alpha / 2
  qnorm(1 - alpha_eff)
}

# -----------------------------------------------------------------------------
# Paired design: discordance / correlation helpers
# -----------------------------------------------------------------------------

# `x` must be a single non-missing number (no range restriction).
check_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x))
    stop(name, " must be a single number.", call. = FALSE)
  invisible(x)
}

# Range of the discordance for marginals p1, p2 (Stark et al., 2022, Eq. 10):
# |p1 - p2| <= psi <= p1 + p2 - 2 p1 p2. The lower end is the smallest possible
# discordance (maximal correlation), the upper end the discordance at rho = 0,
# so the range corresponds to non-negative correlation.
psi_range <- function(p1, p2) {
  c(abs(p1 - p2), p1 + p2 - 2 * p1 * p2)
}

psi_to_rho_raw <- function(p1, p2, psi) {
  p11 <- (p1 + p2 - psi) / 2
  (p11 - p1 * p2) / sqrt(p1 * (1 - p1) * p2 * (1 - p2))
}

# Range of the correlation: 0 (psi maximal) up to the value at minimal psi
rho_range <- function(p1, p2) {
  c(0, psi_to_rho_raw(p1, p2, psi_range(p1, p2)[1]))
}

# Discordance of one endpoint from either psi or rho (at most one of them).
# If neither is given, the lower end of the Stark et al. (2022) range is used.
# `suffix` is "d" or "nd" and only used for argument names in messages.
resolve_psi <- function(psi, rho, p1, p2, suffix) {
  psi_nm <- paste0("psi.", suffix)
  rho_nm <- paste0("rho.", suffix)
  if (!is.null(psi) && !is.null(rho))
    stop("Specify only one of ", psi_nm, " and ", rho_nm, ".", call. = FALSE)
  if (!is.null(rho)) return(rho_to_psi(p1, p2, rho, name = rho_nm))
  if (!is.null(psi)) return(psi)
  psi <- psi_range(p1, p2)[1]
  message(psi_nm, " and ", rho_nm, " not given: using the lower bound of ",
          "Stark et al. (2022), ", psi_nm, " = |p1 - p2| = ", round(psi, 4),
          ". This is the smallest possible discordance and gives the smallest ",
          "sample size. Larger discordance needs more subjects.")
  psi
}
