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
#'     (non-diseased) can be specified, or alternatively the correlations
#'     \code{rho.d} and \code{rho.nd} (one of the two per endpoint). Without
#'     either, the lower Stark bound is used.
#' }
#' Superiority and non-inferiority hypotheses may be specified independently
#' for sensitivity and specificity. The study is successful only if both
#' endpoint null hypotheses are rejected (intersection-union test, IUT), so the
#' type I error rate is controlled at level \code{alpha} without multiplicity
#' adjustment.
#'
#' @details
#' Let \eqn{\pi}{pi} denote the prevalence. Sensitivity is evaluated in the
#' \eqn{N\pi}{N * pi} diseased and specificity in the
#' \eqn{N(1-\pi)}{N * (1 - pi)} non-diseased
#' subjects (per arm in the unpaired design, where \eqn{N}{N} is the size of each
#' arm). These subgroups are independent, so
#' \deqn{\mathrm{Power}_{\mathrm{Overall}} = \mathrm{Power}_{\mathrm{Se}} \times
#'   \mathrm{Power}_{\mathrm{Sp}}.}{Power_Overall = Power_Se * Power_Sp.}
#'
#' \strong{Unpaired design.} Large-sample normal approximation for the
#' difference of two independent proportions with variance
#' \eqn{\theta_C(1-\theta_C) + \theta_E(1-\theta_E)}{theta_C (1 - theta_C) + theta_E (1 - theta_E)}
#' per subject, where \eqn{\theta_C}{theta_C} and \eqn{\theta_E}{theta_E} are
#' the comparator and experimental sensitivity (or specificity), i.e.
#' \code{se.c}, \code{se.e} (or \code{sp.c}, \code{sp.e}). For superiority the
#' critical value uses the null variance chosen by \code{variance}. For
#' non-inferiority the alternative variance is used with the margin
#' \code{delta}. The sample size is the smallest integer \eqn{N}{N} whose overall
#' power reaches the target.
#'
#' \strong{Paired design.} Superiority uses the asymptotic formula of
#' Stark et al. (2022, Miettinen). Non-inferiority uses the restricted maximum
#' likelihood (score) test of Liu et al. (2002), with the variance of the
#' estimated difference evaluated at the boundary of the null hypothesis.
#' Both depend on the discordance \eqn{\psi}{psi}, which must lie in the interval of
#' Stark et al. (2022, Eq. 10),
#' \eqn{|d| \le \psi \le c + e - 2ce}{|d| <= psi <= c + e - 2 * c * e}, where
#' \eqn{c}{c} and \eqn{e}{e} are the comparator and experimental sensitivity
#' (or specificity) and \eqn{d = e - c}{d = e - c} their difference. This
#' corresponds to non-negative correlation between the two tests. Negative
#' correlation is not considered. The sample size increases with \eqn{\psi}{psi}. If
#' no discordance is specified, its lower bound \eqn{|d|}{|d|} is used, as in the
#' worked example of Stark et al. As for the unpaired design, the sample size is
#' the smallest integer \eqn{N}{N} whose overall power reaches the target.
#'
#' The dependence between the tests may be given as a correlation
#' (\code{rho.d}, \code{rho.nd}, phi coefficient of the two binary results)
#' instead of the discordance. The correlation is converted once, at the
#' alternative, to the discordance with \code{\link{rho_to_psi}}. All
#' calculations, including the null quantities of the Liu et al. (2002) score
#' test, then use this discordance exactly as if \code{psi.d} / \code{psi.nd}
#' had been given.
#' The correlation must lie between 0 and its maximum for the given
#' sensitivities (or specificities). Use \code{\link{paired_bounds}} to see the
#' feasible ranges before the calculation.
#'
#' \strong{Exact power (unpaired design).} With \code{exact = TRUE} the power
#' is the exact probability, over two independent binomial counts per
#' endpoint, that the Wald statistic with estimated variance exceeds the
#' critical value (the same test as in \code{\link{power.copr.comp.emp}}).
#' The \code{"null"} variance has no exact counterpart, so \code{"unpooled"}
#' is used unless \code{variance = "pooled"} is given. The exact power is not
#' monotone in \eqn{N}{N}, so the sample size is found by a scan around the
#' asymptotic solution (see \code{exact_rule} and \code{\link{power.copr.test}}).
#' The sums neglect binomial tail probabilities below \eqn{10^{-12}}{1e-12}.
#'
#' @references
#' Stark, M., Hesse, M., Brannath, W., & Zapf, A. (2022).
#' Blinded sample size re-estimation in a comparative diagnostic
#' accuracy study.
#' \emph{BMC Medical Research Methodology}, 22, 115.
#' \doi{10.1186/s12874-022-01564-2}
#'
#' Liu, J., Hsueh, H., Hsieh, E., & Chen, J. J. (2002).
#' Tests for equivalence or non-inferiority for paired binary data.
#' \emph{Statistics in Medicine}, 21(2), 231-245.
#' \doi{10.1002/sim.1012}
#'
#' @param N Sample size, or \code{NULL} to calculate it. Participants per arm
#'   (unpaired) or in total (paired).
#' @param power Desired overall power, or \code{NULL} to calculate it.
#' @param prev Disease prevalence (0 < prev < 1).
#' @param design \code{"unpaired"} (default) or \code{"paired"}.
#' @param se.c Sensitivity of the comparator test.
#' @param se.e Sensitivity of the experimental test.
#' @param sp.c Specificity of the comparator test.
#' @param sp.e Specificity of the experimental test.
#' @param alpha Type I error level per endpoint (default 0.025). It is the
#'   one-sided level unless \code{alternative = "two.sided"}, in which case
#'   the critical value uses \code{alpha / 2}.
#' @param alternative \code{"one.sided"} (default) or \code{"two.sided"}.
#' @param exact Only for \code{design = "unpaired"}: use the exact binomial
#'   power of the Wald test (estimated variance, \code{variance = "unpooled"}
#'   by default or \code{"pooled"}) instead of the normal approximation. \code{N}
#'   must then be an integer, split into \code{round(N * prev)} diseased and
#'   the remaining non-diseased subjects \emph{per arm}. See Details. Default
#'   \code{FALSE}.
#' @param exact_rule Only for \code{exact = TRUE} and a given \code{power}:
#'   \code{"first"} (default, smallest \eqn{N} reaching the target) or
#'   \code{"stays_above"} (smallest \eqn{N} from which the power stays at or
#'   above the target).
#' @param psi.d Discordance probability in the diseased (only for
#'   \code{design = "paired"}). At most one of \code{psi.d} and \code{rho.d} may
#'   be given. If neither is, the lower bound of Stark et al. (2022),
#'   \code{abs(se.c - se.e)}, is used and a message is printed. This is the smallest
#'   possible discordance and gives the smallest sample size.
#' @param psi.nd Discordance probability in the non-diseased. As \code{psi.d},
#'   with default \code{abs(sp.c - sp.e)}.
#' @param rho.d Correlation of the two tests' results in the diseased, as an
#'   alternative to \code{psi.d} (only for \code{design = "paired"}).
#' @param rho.nd Correlation of the two tests' results in the non-diseased, as
#'   an alternative to \code{psi.nd} (only for \code{design = "paired"}).
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
#'   \item{n}{Sample size \code{N}: per arm for the unpaired design, total
#'     number of participants for the paired design.}
#'   \item{n_dis}{Number of diseased subjects (per arm for the unpaired
#'     design).}
#'   \item{n_nondis}{Number of non-diseased subjects (per arm for the unpaired
#'     design).}
#'   \item{prev}{Disease prevalence.}
#'   \item{sig.level}{Type I error level.}
#'   \item{alternative}{\code{"one.sided"} or \code{"two.sided"}.}
#'   \item{power}{Overall power.}
#'   \item{power_se}{Power for the sensitivity endpoint.}
#'   \item{power_sp}{Power for the specificity endpoint.}
#'   \item{type.se}{Hypothesis type for sensitivity.}
#'   \item{type.sp}{Hypothesis type for specificity.}
#'   \item{psi.d, psi.nd}{Discordance used (paired design only). It is converted
#'     from \code{rho.d} / \code{rho.nd} if those were given.}
#'   \item{variance}{Variance used (unpaired design only).}
#'   \item{exact}{Whether exact binomial power was used (unpaired design only).}
#'   \item{target_power}{Requested power (only when \code{power} was given).}
#'   \item{exact_rule}{Rule used (only for \code{exact = TRUE} and a given \code{power}).}
#' }
#'
#' @seealso \code{\link{power.copr.test}} for the single-test design,
#'   \code{\link{paired_bounds}} and \code{\link{rho_to_psi}} for the paired
#'   dependence parameters.
#' @author Maia Muresan & Ainesh Sewak (DCR Bern)
#' @importFrom stats pnorm qnorm qbinom dbinom
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
#' # Unpaired, exact binomial power at an integer N
#' power.copr.comp.test(N = 300, prev = 0.6, exact = TRUE,
#'                      se.c = 0.75, se.e = 0.85,
#'                      sp.c = 0.65, sp.e = 0.75)
#'
#' # Unpaired, specificity non-inferiority
#' power.copr.comp.test(power = 0.80, prev = 0.6,
#'                      se.c = 0.75, se.e = 0.85,
#'                      sp.c = 0.65, sp.e = 0.75,
#'                      delta.sp = 0.05, type.sp = "ni")
#'
#' # Paired design with discordances (default if omitted: lower Stark bound)
#' power.copr.comp.test(power = 0.80, prev = 0.47, design = "paired",
#'                      se.c = 0.81, se.e = 0.90,
#'                      sp.c = 0.66, sp.e = 0.80,
#'                      psi.d = 0.09, psi.nd = 0.14)
#' power.copr.comp.test(N = 134, prev = 0.47, design = "paired",
#'                      se.c = 0.81, se.e = 0.90,
#'                      sp.c = 0.66, sp.e = 0.80,
#'                      psi.d = 0.09, psi.nd = 0.14,
#'                      delta.sp = 0.05, type.sp = "ni")
#'
#' # Paired design with correlations instead of discordances
#' paired_bounds(se.c = 0.81, se.e = 0.90, sp.c = 0.66, sp.e = 0.80)
#' power.copr.comp.test(power = 0.80, prev = 0.47, design = "paired",
#'                      se.c = 0.81, se.e = 0.90,
#'                      sp.c = 0.66, sp.e = 0.80,
#'                      rho.d = 0.5, rho.nd = 0.3)
power.copr.comp.test <- function(N = NULL,
                                 power = NULL,
                                 prev,
                                 design = c("unpaired", "paired"),
                                 se.c, se.e,
                                 sp.c, sp.e,
                                 alpha = 0.025,
                                 psi.d = NULL,
                                 psi.nd = NULL,
                                 rho.d = NULL,
                                 rho.nd = NULL,
                                 delta.se = 0,
                                 delta.sp = 0,
                                 type.se = c("sup", "ni"),
                                 type.sp = c("sup", "ni"),
                                 variance = c("null", "pooled", "unpooled"),
                                 alternative = c("one.sided", "two.sided"),
                                 exact = FALSE,
                                 exact_rule = c("first", "stays_above")) {

  design      <- match.arg(design)
  type.se     <- match.arg(type.se)
  type.sp     <- match.arg(type.sp)
  alternative <- match.arg(alternative)
  exact_rule  <- match.arg(exact_rule)
  if (!missing(variance) && design == "paired")
    warning("'variance' is only used for design = \"unpaired\" and is ignored.")
  var_missing <- missing(variance)
  variance <- match.arg(variance)

  # -----------------------------
  # Input validation
  # -----------------------------
  check_one_of(N, power)
  check_inputs(N = N, power = power, prev = prev, alpha = alpha,
               probs = list(se.c = se.c, se.e = se.e, sp.c = sp.c, sp.e = sp.e))
  warn_alpha(alpha, alternative)
  check_num(delta.se, "delta.se", -1, 1)
  check_num(delta.sp, "delta.sp", -1, 1)

  if (exact) {
    if (design == "paired")
      stop("exact = TRUE is only available for design = \"unpaired\".")
    if (!is.null(N) && abs(N - round(N)) > 1e-8)
      stop("N must be an integer when exact = TRUE.")
    # the exact test estimates the variance from the data, so the "null"
    # variance (true comparator proportion) has no exact counterpart
    if (var_missing) variance <- "unpooled"
    if (variance == "null")
      stop("exact = TRUE requires variance = \"pooled\" or \"unpooled\".")
  }

  # one-sided level used for the critical value
  alpha_eff <- if (alternative == "one.sided") alpha else alpha / 2

  if (design == "paired") {
    # the correlation is converted once, at the alternative; everything below
    # works with the discordance (default: lower Stark bound)
    psi.d  <- resolve_psi(psi.d,  rho.d,  se.c, se.e, "d")
    psi.nd <- resolve_psi(psi.nd, rho.nd, sp.c, sp.e, "nd")
    check_paired(se.c, se.e, sp.c, sp.e, psi.d, psi.nd,
                 delta.se, delta.sp, type.se, type.sp)
  } else {
    if (!is.null(psi.d) || !is.null(psi.nd) || !is.null(rho.d) || !is.null(rho.nd))
      warning("'psi.d', 'psi.nd', 'rho.d' and 'rho.nd' are only used for ",
              "design = \"paired\" and are ignored.")
    check_unpaired(se.c, se.e, sp.c, sp.e, delta.se, delta.sp,
                   type.se, type.sp)
  }

  # power of one endpoint with n subjects in its subgroup
  endpoint_power <- function(n, p.c, p.e, delta, type, psi, use_exact) {
    if (use_exact) {
      unpaired_exact_endpoint_power(n, alpha_eff, p.c, p.e, delta, type, variance)
    } else if (design == "unpaired") {
      unpaired_endpoint_power(n, alpha_eff, p.c, p.e, delta, type, variance)
    } else {
      paired_endpoint_power(n, alpha_eff, p.c, p.e, psi, delta, type)
    }
  }

  # subgroup sizes at total N; exact uses round(N * prev) diseased and the
  # remaining non-diseased
  n_diseased <- function(N, use_exact) if (use_exact) round(N * prev) else N * prev

  # endpoint powers (se, sp) at total sample size N (asymptotic unless
  # use_exact; the asymptotic version is also the starting point of exact_n())
  endpoint_powers <- function(N, use_exact = FALSE) {
    n_dis <- n_diseased(N, use_exact)
    c(se = endpoint_power(n_dis,     se.c, se.e, delta.se, type.se, psi.d, use_exact),
      sp = endpoint_power(N - n_dis, sp.c, sp.e, delta.sp, type.sp, psi.nd, use_exact))
  }

  # -----------------------------
  # Sample size
  # -----------------------------
  if (is.null(N)) {
    N <- find_n(function(N) prod(endpoint_powers(N)), power)
    if (exact) {
      # exact power is a sawtooth in N: scan around the asymptotic solution
      exact_pw <- function(N) vapply(N, function(n) prod(endpoint_powers(n, TRUE)),
                                     numeric(1))
      N <- exact_n(exact_pw, power, N_asym = N, rule = exact_rule)
    }
  }

  pw <- endpoint_powers(N, exact)
  n_dis <- n_diseased(N, exact)

  # -----------------------------
  # Return structured object
  # -----------------------------
  if (design == "unpaired") {
    method_design <- "Comparative unpaired design"
    approx <- if (exact) {
      paste0("Exact binomial (Wald test, ", variance, " variance)")
    } else {
      "Wald normal approximation"
    }
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
    n = N, n_dis = n_dis, n_nondis = N - n_dis, prev = prev,
    sig.level = alpha, alternative = alternative,
    power = prod(pw), power_se = pw[["se"]], power_sp = pw[["sp"]],
    target_power = power,
    type.se = type.se, type.sp = type.sp,
    psi.d = if (design == "paired") psi.d,
    psi.nd = if (design == "paired") psi.nd,
    variance = if (design == "unpaired") variance,
    exact = if (design == "unpaired") exact,
    exact_rule = if (exact && !is.null(power)) exact_rule
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

  # Discordance range of Stark et al. (2022, Eq. 10): |d| <= psi <= c + e - 2ce
  psi_bounds <- function(psi, c, e, name) {
    lower <- psi_range(c, e)[1]
    upper <- psi_range(c, e)[2]
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
  p.bar <- (c + e) / 2
  V0 <- switch(variance,
    pooled   = 2 * p.bar * (1 - p.bar),
    unpooled = VA,
    null     = if (type == "sup") v(c, c) else VA)

  1 - pnorm(z.alpha * sqrt(V0 / n), (e - c) + margin, sqrt(VA / n))
}

# Score quantities for paired non-inferiority (Liu et al., 2002, Stat Med 21:
# 231-245, Section 3.2): w_u (the paper's w_l, ratio of the true and the
# restricted-MLE standard error) and sigma^2 * n.
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

# Exact power of one endpoint, unpaired design, for n per arm (integer):
# P(reject) under X ~ Bin(n, c), Y ~ Bin(n, e) independent, where the Wald
# statistic (y/n - x/n + margin) / se_hat > z uses the unpooled or pooled
# variance estimate (as in power.copr.comp.emp()). Zero/zero counts as not
# rejected. The sums run over the 1e-12 .. 1 - 1e-12 binomial quantile ranges,
# so the neglected probability is below 4e-12.
unpaired_exact_endpoint_power <- function(n, alpha, c, e, delta, type, variance) {
  if (n < 1) return(0)
  z <- qnorm(1 - alpha)
  margin <- if (type == "ni") delta else 0
  rng <- function(p) qbinom(1e-12, n, p):qbinom(1 - 1e-12, n, p)
  x <- rng(c)
  y <- rng(e)
  px <- dbinom(x, n, c)
  py <- dbinom(y, n, e)
  xh <- x / n
  yh <- y / n

  # rows: y, columns: x
  num <- outer(yh, xh, "-") + margin
  se <- if (variance == "pooled") {
    pb <- outer(yh, xh, "+") / 2
    sqrt(2 * pb * (1 - pb) / n)
  } else {
    sqrt(outer(yh * (1 - yh), xh * (1 - xh), "+") / n)
  }
  rej <- num / se > z
  rej[is.na(rej)] <- FALSE
  sum(outer(py, px) * rej)
}
