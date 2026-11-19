#' Power calculation for a single-test diagnostic accuracy study with co-primary endpoints
#'
#' @description
#' The function determines the achieved overall power for a confirmatory diagnostic
#' accuracy study with co-primary endpoints (sensitivity and specificity)
#' for a fixed total sample size.
#'
#' Both endpoints are tested for superiority against fixed null values.
#' Study success requires rejection of both null hypotheses.
#' Overall power is therefore defined under an intersection–union
#' testing (IUT) framework, meaning that both tests must be significant.
#' Under this framework, the global Type I error rate is controlled
#' at level \eqn{\alpha}{alpha} without multiplicity adjustment.
#'
#' The calculation explicitly incorporates the assumed disease prevalence,
#' which determines the expected number of diseased and non-diseased subjects.
#'
#' @details
#' Let \eqn{\pi}{pi} denote the assumed disease prevalence. The study population
#' is partitioned into diseased (\eqn{N_{Se} = N\pi}{N_Se = N * pi}) and non-diseased
#' (\eqn{N_{Sp} = N(1-\pi)}{N_Sp = N * (1 - pi)}) independent subgroups.
#'
#' Endpoint-specific powers are derived from the asymptotic normal
#' test statistics corresponding to the sample size formulas.
#' The implementation uses large-sample approximations for binomial
#' proportions based on either:
#' \itemize{
#'   \item the Wald approximation on the probability scale, or
#'   \item the logit-transformed normal approximation.
#' }
#'
#' Under standard binomial sampling assumptions, the estimators of
#' sensitivity and specificity are asymptotically independent.
#'
#' Overall power is the probability that both test statistics exceed their
#' critical values. Because the endpoints are independent, this is
#'
#' \deqn{\mathrm{Power}_{\mathrm{overall}} = \mathrm{Power}_{\mathrm{Se}} \times
#'   \mathrm{Power}_{\mathrm{Sp}}.}{Power_overall = Power_Se * Power_Sp.}
#'
#' The function returns the achieved endpoint-specific powers and
#' their joint probability for the specified total sample size.
#' The power calculation is consistent with the asymptotic normal framework
#' underlying the sample size methodology of Stark and Zapf (2020).
#' @references
#' Stark, M., & Zapf, A. (2020). Sample size calculation and re-estimation
#' based on the prevalence in a single-arm confirmatory diagnostic accuracy
#' study. \emph{Statistical Methods in Medical Research}, 29(10), 2958–2971.
#' \doi{10.1177/0962280220913588}
#'
#' @param N Total sample size.
#' @param prev Disease prevalence (0 < prev < 1).
#' @param se.e Sensitivity under the alternative hypothesis.
#' @param se.min Sensitivity under the null hypothesis.
#' @param sp.e Specificity under the alternative hypothesis.
#' @param sp.min Specificity under the null hypothesis.
#' @param alpha Type I error level.
#' @param ci_method Normal approximation method: "wald" (default) or "logit".
#' @param alternative "one.sided" (default) or "two.sided".
#' @param exact If \code{TRUE}, the endpoint powers are the exact binomial
#'   rejection probabilities of the Wald/logit tests instead of the normal
#'   approximation. \code{N} must then be an integer, split into
#'   \code{round(N * prev)} diseased and the remaining non-diseased
#'   individuals. Default \code{FALSE}.
#' @author Maia Muresan (DCR Bern)
#'
#' @importFrom stats pnorm pbinom plogis qlogis
#' @return An object of class \code{"power.htest"}, a list with the elements:
#' \describe{
#'   \item{n}{Total sample size.}
#'   \item{n_dis}{Total number of diseased subjects.}
#'   \item{n_nondis}{Total number of non-diseased subjects.}
#'   \item{prev}{Disease prevalence.}
#'   \item{sig.level}{Type I error level.}
#'   \item{alternative}{\code{"one.sided"} or \code{"two.sided"}.}
#'   \item{power}{Overall power.}
#'   \item{power_se}{Power for the sensitivity endpoint.}
#'   \item{power_sp}{Power for the specificity endpoint.}
#'   \item{ci_method}{Normal approximation used.}
#'   \item{exact}{Whether exact binomial power was used.}
#' }
#'
#' @noRd
#' @examples
#'
#' # Example 1: Achieved power for a fixed sample size
#' # Objective: evaluate power to demonstrate that sensitivity
#' # and specificity exceed pre-specified minimum values
#' power_opt_bin(
#'   N = 122,            # total sample size
#'   prev = 0.5,         # disease prevalence
#'   se.e = 0.85,        # expected sensitivity
#'   se.min = 0.65,      # minimum acceptable sensitivity
#'   sp.e = 0.85,        # expected specificity
#'   sp.min = 0.55,      # minimum acceptable specificity
#'   alpha = 0.025       # one-sided significance level
#' )
#' #same
#' power_opt_bin(
#'   N = 122,            # total sample size
#'   prev = 0.5,         # disease prevalence
#'   se.e = 0.85,        # expected sensitivity
#'   se.min = 0.65,      # minimum acceptable sensitivity
#'   sp.e = 0.85,        # expected specificity
#'   sp.min = 0.55,      # minimum acceptable specificity
#'   alpha = 0.025,      # one-sided significance level (default)
#'
#' )
#'
#' # Example 2: Same design using logit-based confidence interval
#' power_opt_bin(
#'   N = 122,
#'   prev = 0.5,
#'   se.e = 0.85,
#'   se.min = 0.65,
#'   sp.e = 0.85,
#'   sp.min = 0.55,
#'   alpha = 0.025,
#'   ci_method = "logit",
#'   alternative = "one.sided"
#' )

power_opt_bin <- function(N,
                          prev,
                          se.e, se.min,
                          sp.e, sp.min,
                          alpha = 0.025, ci_method=c("wald", "logit"),
                          alternative=c("one.sided","two.sided"),
                          exact = FALSE
) {

  ci_method   <- match.arg(ci_method)
  alternative <- match.arg(alternative)
  # -----------------------------
  # Input validation
  # -----------------------------
  check_inputs(N = N, prev = prev, alpha = alpha)
  check_single(se.e, se.min, sp.e, sp.min)
  if (exact && abs(N - round(N)) > 1e-8)
    stop("N must be an integer when exact = TRUE.")

  # -----------------------------
  # Common quantities
  # -----------------------------
  #critical value
  z_alpha <- z_critical(alpha, alternative)

  if (exact) {
    # exact binomial power with integer subgroup sizes
    ex <- exact_power_bin(round(N), prev, se.e, se.min, sp.e, sp.min,
                          z_alpha, ci_method)
    Nd <- ex$n_dis
    Nh <- ex$n_nondis
    power_se <- ex$power_se
    power_sp <- ex$power_sp
  } else {
    #sample sizes per subgroup
    Nd <- N * prev           # number of diseased
    Nh <- N * (1 - prev)     # number of non-diseased

    # -----------------------------
    # Scale-specific quantities (probability or logit scale)
    # -----------------------------
    sc_se <- endpoint_scale(se.e, se.min, ci_method)
    sc_sp <- endpoint_scale(sp.e, sp.min, ci_method)
    d_se     <- sc_se$d
    v_se.min <- sc_se$v_min
    v_se.e   <- sc_se$v_e
    d_sp     <- sc_sp$d
    v_sp.min <- sc_sp$v_min
    v_sp.e   <- sc_sp$v_e

    # -----------------------------
    # Endpoint-specific power
    # -----------------------------
    #Sensitivity power
    z_beta_se <- (sqrt(Nd) * d_se - z_alpha * sqrt(v_se.min)) / sqrt(v_se.e)
    power_se <- pnorm(z_beta_se)

    # Specificity power
    z_beta_sp <- (sqrt(Nh) * d_sp - z_alpha * sqrt(v_sp.min)) / sqrt(v_sp.e)
    power_sp <- pnorm(z_beta_sp)
  }

  # -----------------------------
  # Overall power (IUT framework)
  # -----------------------------
  # sensitivity and specificity come from disjoint samples (independent)
  overall_power <- power_se * power_sp


  # -----------------------------
  # Return structured object
  # -----------------------------
  power_htest(
    method = power_method("Single-test design",
                          single_approx(ci_method, exact), "power"),
    n = N, n_dis = Nd, n_nondis = Nh, prev = prev,
    sig.level = alpha, alternative = alternative,
    power = overall_power, power_se = power_se, power_sp = power_sp,
    ci_method = ci_method, exact = exact
  )
}


