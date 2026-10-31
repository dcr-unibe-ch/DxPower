#' Power calculations for co-primary diagnostic accuracy endpoints - Comparative designs
#'
#' @description
#' Computes either the required total sample size or the achieved overall power
#' of a comparative confirmatory diagnostic accuracy study with co-primary
#' endpoints (sensitivity and specificity). Exactly one of \code{N} and
#' \code{power} must be \code{NULL}. The one that is \code{NULL} is calculated.
#'
#' Two designs are available:
#' \itemize{
#'   \item \code{design = "unpaired"}: subjects are randomised to the
#'     experimental or the comparator test. The endpoints compare two
#'     independent proportions.
#'   \item \code{design = "paired"}: all subjects receive both tests. The
#'     discordance probabilities \code{psi.d} (diseased) and \code{psi.nd}
#'     (non-diseased) must be specified.
#' }
#' Superiority and non-inferiority hypotheses may be specified independently
#' for sensitivity and specificity. The study is successful only if both
#' endpoint null hypotheses are rejected (intersection-union test, IUT), so the
#' type I error rate is controlled at level \code{alpha} without multiplicity
#' adjustment.
#'
#' @details
#' Let \eqn{\pi} denote the prevalence. Sensitivity is evaluated in the
#' \eqn{N\pi} diseased and specificity in the \eqn{N(1-\pi)} non-diseased
#' subjects. These subgroups are independent, so
#' \deqn{Power_{overall} = Power_{Se} \times Power_{Sp}.}
#'
#' \strong{Unpaired design.} Large-sample normal approximation for the
#' difference of two independent proportions with variance
#' \eqn{\theta_C(1-\theta_C) + \theta_E(1-\theta_E)}. For superiority the
#' critical value uses the null variance chosen by \code{variance}. For
#' non-inferiority the alternative variance is used with the margin
#' \code{delta}. The sample size is the smallest integer \eqn{N} whose overall
#' power reaches the target.
#'
#' \strong{Paired design.} Superiority uses the asymptotic formula of
#' Stark et al. (2022, Miettinen). Non-inferiority uses the Liu score method.
#' Both depend on the discordance \eqn{\psi}, which must satisfy
#' \eqn{|d| \le \psi \le c + e - 2ce}. The sample size splits the allowed type
#' II error between the endpoints numerically such that both endpoints imply
#' the same total sample size and the overall power equals the target. The
#' returned \code{n} is not rounded.
#'
#' @references
#' Stark, M., Hesse, M., Brannath, W., & Zapf, A. (2022).
#' Blinded sample size re-estimation in a comparative diagnostic
#' accuracy study.
#' \emph{BMC Medical Research Methodology}, 22, 115.
#' \doi{10.1186/s12874-022-01564-2}
#'
#' @param N Total sample size, or \code{NULL} to calculate it.
#' @param power Desired overall power, or \code{NULL} to calculate it.
#' @param prev Disease prevalence (0 < prev < 1).
#' @param design \code{"unpaired"} (default) or \code{"paired"}.
#' @param se.c Sensitivity of the comparator test.
#' @param se.e Sensitivity of the experimental test.
#' @param sp.c Specificity of the comparator test.
#' @param sp.e Specificity of the experimental test.
#' @param alpha One-sided type I error rate per endpoint (default 0.025).
#' @param psi.d Discordance probability in the diseased. Required (and only
#'   used) for \code{design = "paired"}.
#' @param psi.nd Discordance probability in the non-diseased. Required (and
#'   only used) for \code{design = "paired"}.
#' @param delta.se Non-inferiority margin for sensitivity
#'   (used only if \code{type.se = "ni"}).
#' @param delta.sp Non-inferiority margin for specificity
#'   (used only if \code{type.sp = "ni"}).
#' @param type.se Type of hypothesis for sensitivity:
#'   \code{"sup"} (superiority, default) or \code{"ni"} (non-inferiority).
#' @param type.sp Type of hypothesis for specificity:
#'   \code{"sup"} (superiority, default) or \code{"ni"} (non-inferiority).
#' @param variance Only for \code{design = "unpaired"}: variance of the
#'   difference estimator under the null used for the critical value:
#'   \code{"null"} (default, variance at the comparator proportion for
#'   superiority, alternative variance for non-inferiority), \code{"pooled"}
#'   (pooled proportion \eqn{(c+e)/2}) or \code{"unpooled"} (alternative
#'   variance). The alternative variance is always unpooled.
#'
#' @return An object of class \code{"power.htest"}, a list with the elements:
#' \describe{
#'   \item{n}{Total sample size.}
#'   \item{n_dis}{Total number of diseased subjects.}
#'   \item{n_nondis}{Total number of non-diseased subjects.}
#'   \item{prev}{Disease prevalence.}
#'   \item{sig.level}{Type I error level.}
#'   \item{alternative}{\code{"one.sided"}.}
#'   \item{power}{Overall power.}
#'   \item{power_se}{Power for the sensitivity endpoint.}
#'   \item{power_sp}{Power for the specificity endpoint.}
#'   \item{type.se}{Hypothesis type for sensitivity.}
#'   \item{type.sp}{Hypothesis type for specificity.}
#'   \item{target_power}{Requested power (only when \code{power} was given).}
#' }
#'
#' @seealso \code{\link{power.copr.test}} for the single-test design.
#' @author Maia Muresan & Ainesh Sewak (DCR Bern)
#' @importFrom stats pnorm qnorm
#' @export
#' @examples
#' # Unpaired design: sample size, then power at that sample size
#' power.copr.comp.test(power = 0.80, prev = 0.6,
#'                      se.c = 0.75, se.e = 0.85,
#'                      sp.c = 0.65, sp.e = 0.75)
#' power.copr.comp.test(N = 910, prev = 0.6,
#'                      se.c = 0.75, se.e = 0.85,
#'                      sp.c = 0.65, sp.e = 0.75)
#'
#' # Unpaired, specificity non-inferiority
#' power.copr.comp.test(power = 0.80, prev = 0.6,
#'                      se.c = 0.75, se.e = 0.85,
#'                      sp.c = 0.65, sp.e = 0.75,
#'                      delta.sp = 0.05, type.sp = "ni")
#'
#' # Paired design: psi.d and psi.nd are required
#' power.copr.comp.test(power = 0.80, prev = 0.47, design = "paired",
#'                      se.c = 0.81, se.e = 0.90,
#'                      sp.c = 0.66, sp.e = 0.80,
#'                      psi.d = 0.09, psi.nd = 0.14)
#' power.copr.comp.test(N = 134, prev = 0.47, design = "paired",
#'                      se.c = 0.81, se.e = 0.90,
#'                      sp.c = 0.66, sp.e = 0.80,
#'                      psi.d = 0.09, psi.nd = 0.14,
#'                      delta.sp = 0.05, type.sp = "ni")
power.copr.comp.test <- function(N = NULL,
                                 power = NULL,
                                 prev,
                                 design = c("unpaired", "paired"),
                                 se.c, se.e,
                                 sp.c, sp.e,
                                 alpha = 0.025,
                                 psi.d = NULL,
                                 psi.nd = NULL,
                                 delta.se = 0,
                                 delta.sp = 0,
                                 type.se = c("sup", "ni"),
                                 type.sp = c("sup", "ni"),
                                 variance = c("null", "pooled", "unpooled")) {

  design   <- match.arg(design)
  type.se  <- match.arg(type.se)
  type.sp  <- match.arg(type.sp)
  if (!missing(variance) && design == "paired")
    warning("'variance' is only used for design = \"unpaired\" and is ignored.")
  variance <- match.arg(variance)

  # -----------------------------
  # Input validation
  # -----------------------------
  if (sum(vapply(list(N, power), is.null, logical(1))) != 1) {
    stop("exactly one of 'N' and 'power' must be NULL")
  }
  if (!is.null(N) && N <= 0) stop("N must be positive.")
  if (!is.null(power) && (power <= 0 || power >= 1))
    stop("power must be between 0 and 1.")
  if (prev <= 0 || prev >= 1) stop("Prevalence must be between 0 and 1.")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1.")

  if (design == "paired") {
    if (is.null(psi.d) || is.null(psi.nd))
      stop("psi.d and psi.nd must be specified for design = \"paired\".")
    check_paired(se.c, se.e, sp.c, sp.e, psi.d, psi.nd,
                 delta.se, delta.sp, type.se, type.sp)
  } else {
    if (!is.null(psi.d) || !is.null(psi.nd))
      warning("'psi.d' and 'psi.nd' are only used for design = \"paired\" ",
              "and are ignored.")
    check_unpaired(se.c, se.e, sp.c, sp.e, delta.se, delta.sp,
                   type.se, type.sp)
  }

  # endpoint powers (se, sp) at total sample size N
  endpoint_powers <- function(N) {
    n_dis    <- N * prev
    n_nondis <- N * (1 - prev)
    if (design == "unpaired") {
      c(se = unpaired_endpoint_power(n_dis, alpha, se.c, se.e, delta.se,
                                     type.se, variance),
        sp = unpaired_endpoint_power(n_nondis, alpha, sp.c, sp.e, delta.sp,
                                     type.sp, variance))
    } else {
      c(se = paired_endpoint_power(n_dis, alpha, se.c, se.e, psi.d, delta.se,
                                   type.se),
        sp = paired_endpoint_power(n_nondis, alpha, sp.c, sp.e, psi.nd,
                                   delta.sp, type.sp))
    }
  }

  # -----------------------------
  # Sample size
  # -----------------------------
  if (is.null(N)) {
    N <- if (design == "unpaired") {
      find_n(function(N) prod(endpoint_powers(N)), power)
    } else {
      paired_n(power, prev, alpha, se.c, se.e, sp.c, sp.e, psi.d, psi.nd,
               delta.se, delta.sp, type.se, type.sp)
    }
  }

  pw <- endpoint_powers(N)

  # -----------------------------
  # Return structured object
  # -----------------------------
  if (design == "unpaired") {
    method_design <- "Comparative unpaired design"
    approx <- "Wald normal approximation"
  } else {
    scores <- c(sup = "Miettinen", ni = "Liu score")
    method_design <- "Comparative paired design"
    approx <- paste0("Score normal approximation (",
                     paste(scores[unique(c(type.se, type.sp))],
                           collapse = " / "), ")")
  }

  power_htest(
    method = power_method(method_design, approx,
                          if (is.null(power)) "power" else "samplesize"),
    n = N, n_dis = N * prev, n_nondis = N * (1 - prev), prev = prev,
    sig.level = alpha, alternative = "one.sided",
    power = prod(pw), power_se = pw[["se"]], power_sp = pw[["sp"]],
    target_power = power,
    type.se = type.se, type.sp = type.sp
  )
}

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

# Hypothesis consistency checks shared by both designs.
check_hypotheses <- function(se.c, se.e, sp.c, sp.e, delta.se, delta.sp,
                             type.se, type.sp) {
  # Superiority requires expected improvement
  if (type.se == "sup" && se.e <= se.c)
    stop("Sensitivity superiority requires se.e > se.c")
  if (type.sp == "sup" && sp.e <= sp.c)
    stop("Specificity superiority requires sp.e > sp.c.")
  # Non-inferiority requires true effect above NI boundary
  if (type.se == "ni" && (se.e - se.c) <= -delta.se)
    stop("Sensitivity NI invalid: true difference (se.e - se.c) must be > -delta.se")
  if (type.sp == "ni" && (sp.e - sp.c) <= -delta.sp)
    stop("Specificity NI invalid: true difference (sp.e - sp.c) must be > -delta.sp")
}

check_unpaired <- function(se.c, se.e, sp.c, sp.e, delta.se, delta.sp,
                           type.se, type.sp) {
  check_hypotheses(se.c, se.e, sp.c, sp.e, delta.se, delta.sp, type.se, type.sp)
}

check_paired <- function(se.c, se.e, sp.c, sp.e, psi.d, psi.nd,
                         delta.se, delta.sp, type.se, type.sp) {
  check_hypotheses(se.c, se.e, sp.c, sp.e, delta.se, delta.sp, type.se, type.sp)

  # Informative warning: large NI margin makes the test easy to pass
  if (type.se == "ni" && delta.se > abs(se.e - se.c)) {
    warning("Sensitivity NI margin (delta.se = ", delta.se,
            ") is larger than the expected difference (se.e - se.c = ",
            round(se.e - se.c, 3), "). ",
            "Non-inferiority may be easy to demonstrate with this margin.")
  }
  if (type.sp == "ni" && delta.sp > abs(sp.e - sp.c)) {
    warning("Specificity NI margin (delta.sp = ", delta.sp,
            ") is larger than the expected difference (sp.e - sp.c = ",
            round(sp.e - sp.c, 3), "). ",
            "Non-inferiority may be easy to demonstrate with this margin.")
  }

  # Discordance bounds for paired binary data: |d| <= psi <= c + e - 2ce
  psi_bounds <- function(psi, c, e, name) {
    lower <- abs(c - e)
    upper <- c + e - 2 * c * e
    if (psi < lower || psi > upper) {
      stop(name, " outside theoretical bounds for paired binary data. ",
           name, " must satisfy: ", round(lower, 4), " <= ", name, " <= ",
           round(upper, 4), ". You provided ", name, " = ", psi)
    }
  }
  psi_bounds(psi.d,  se.c, se.e, "psi.d")
  psi_bounds(psi.nd, sp.c, sp.e, "psi.nd")
}

# Power of one endpoint, unpaired design: test "estimated difference > threshold"
#   mean      = true difference + NI margin (0 for superiority)
#   var       = variance of the difference estimator under the alternative / n
#   threshold = z_alpha * sd of the estimator under the null / sqrt(n)
unpaired_endpoint_power <- function(n, alpha, c, e, delta, type, variance) {
  v <- function(c, e) c * (1 - c) + e * (1 - e)
  z.alpha <- qnorm(1 - alpha)
  margin  <- if (type == "ni") delta else 0

  VA <- v(c, e)
  V0 <- switch(variance,
    pooled   = {p.bar <- (c + e) / 2; 2 * p.bar * (1 - p.bar)},
    unpooled = VA,
    null     = if (type == "sup") v(c, c) else VA)

  1 - pnorm(z.alpha * sqrt(V0 / n), (e - c) + margin, sqrt(VA / n))
}

# Liu score quantities for paired non-inferiority: w_u and sigma^2 * n
paired_ni_terms <- function(c, e, psi, delta) {
  d <- e - c
  if (delta + d <= 0)
    stop("NI alternative not satisfied: delta + (e - c) must be > 0.")

  p01 <- (psi - d) / 2
  b_u <- delta * (delta + 1) * p01
  a_u <- (-d) * (1 - delta) - 2 * (p01 + delta)

  disc <- a_u^2 - 8 * b_u
  if (disc <= 0)
    stop("Invalid NI configuration: no real solution exists for the assumed ",
         "margin (delta = ", delta, ") and discordance (psi = ", psi, "). ",
         "The margin may be too small or psi too close to its boundary.")

  p_bar_01 <- (-a_u + sqrt(disc)) / 4
  list(w_u = sqrt(2 * p01 + d - d^2) / sqrt(2 * p_bar_01 - delta - delta^2),
       var = psi - d^2)
}

# Power of one endpoint, paired design
# (sup: Miettinen / Stark 2022, ni: Liu score)
paired_endpoint_power <- function(n, alpha, c, e, psi, delta = 0, type) {
  z.alpha <- qnorm(1 - alpha)
  d <- e - c

  if (type == "sup") {
    z <- (sqrt(n * psi) * d - z.alpha * psi) /
      sqrt(psi^2 - 0.25 * d^2 * (3 + psi))
    pnorm(z)
  } else {
    ni <- paired_ni_terms(c, e, psi, delta)
    z <- (z.alpha / ni$w_u) - (delta + d) / sqrt(ni$var / n)
    1 - pnorm(z)
  }
}

# Required number of subjects of one endpoint, paired design
paired_endpoint_n <- function(alpha, beta, c, e, psi, delta = 0, type) {
  z.beta  <- qnorm(1 - beta)
  z.alpha <- qnorm(1 - alpha)
  d <- e - c

  if (type == "sup") {
    (z.alpha * psi + z.beta * sqrt(psi^2 - 0.25 * d^2 * (3 + psi)))^2 /
      (psi * d^2)
  } else {
    ni <- paired_ni_terms(c, e, psi, delta)
    (ni$var * (z.beta + z.alpha / ni$w_u)^2) / (delta + d)^2
  }
}

# Total sample size of the paired design reaching the overall power `power`.
# The type II error is allocated under the IUT, (1 - beta.se)(1 - beta.sp) =
# power, such that both endpoints imply the same total N.
paired_n <- function(power, prev, alpha, se.c, se.e, sp.c, sp.e, psi.d, psi.nd,
                     delta.se, delta.sp, type.se, type.sp) {
  N_for <- function(beta.se) {
    beta.sp <- 1 - power / (1 - beta.se)
    n.se <- paired_endpoint_n(alpha, beta.se, se.c, se.e, psi.d, delta.se, type.se)
    n.sp <- paired_endpoint_n(alpha, beta.sp, sp.c, sp.e, psi.nd, delta.sp, type.sp)
    c(se = n.se / prev, sp = n.sp / (1 - prev))
  }

  beta.grid <- seq(1e-6, 1 - power - 1e-6, length.out = 10000)
  diff.values <- vapply(beta.grid, function(beta.se) {
    beta.sp <- 1 - power / (1 - beta.se)
    if (beta.sp <= 0 || beta.sp >= 1) return(Inf)
    abs(diff(N_for(beta.se)))
  }, numeric(1))

  max(N_for(beta.grid[which.min(diff.values)]))
}
