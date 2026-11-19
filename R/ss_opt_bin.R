#' Optimal sample size calculation for co-primary diagnostic accuracy endpoints - Single-test design
#'
#' @description
#' The function determines the minimal total sample size required to achieve
#' a prespecified overall power in a confirmatory diagnostic accuracy study with
#' co-primary endpoints (sensitivity and specificity).
#'
#' Both endpoints must demonstrate superiority over fixed null values.
#' Overall power is defined under an intersection–union testing (IUT)
#' framework, meaning that both null hypotheses must be rejected for
#' the study to be successful. Under this framework, the global Type I
#' error rate is controlled at level \eqn{\alpha}{alpha} without multiplicity
#' adjustment.
#'
#' The calculation explicitly accounts for disease prevalence, which
#' determines the expected number of diseased and non-diseased subjects.
#'
#' @details
#' Let \eqn{\pi}{pi} denote the assumed disease prevalence. The study population
#' is partitioned into diseased (\eqn{N_{\mathrm{Se}} = N\pi}{N_Se = N * pi}) and non-diseased
#' (\eqn{N_{\mathrm{Sp}} = N(1-\pi)}{N_Sp = N * (1 - pi)}) independent subgroups.
#'
#' When prevalence is very low or very high, one subgroup may require
#' substantially more subjects than the other. A naive calculation (e.g. taking
#' the maximum of two independently calculated sample sizes) may
#' therefore overpower one endpoint while the other determines the
#' required total sample size.
#'
#' This function implements a joint sample size calculation that accounts
#' for both endpoints simultaneously: the required total sample size is the
#' smallest integer \eqn{N}{N} for which the overall power, as computed by
#' \code{power_opt_bin()}, reaches the target. Neither sensitivity nor
#' specificity is treated as the sole driver of the sample size.
#' Power is based on large-sample normal approximations
#' for binomial proportions using either:
#' \itemize{
#'   \item the Wald approximation on the probability scale, or
#'   \item the logit-transformed normal approximation.
#' }
#'
#' Overall power is the probability that both tests reject. Because the
#' endpoints are independent, this equals
#' \deqn{\mathrm{Power}_{\mathrm{Overall}} = \mathrm{Power}_{\mathrm{Se}} \times
#'   \mathrm{Power}_{\mathrm{Sp}}.}{Power_Overall = Power_Se * Power_Sp.}
#' Without \code{exact}, \eqn{N}{N} is found by doubling an upper bound and
#' then bisection, since the asymptotic power is non-decreasing in \eqn{N}{N}.
#'
#' With \code{exact = TRUE} the power is the exact binomial probability that
#' both tests reject. It is not monotone in \eqn{N}{N} (it is a sawtooth), so the
#' bisection is replaced by a vectorised scan over integer \eqn{N}{N} in a window
#' around the asymptotic solution (0.5 to 1.5 times, extended
#' automatically if needed). Subgroup sizes are \code{round(N * prev)} and the
#' remainder. \code{exact_rule} selects the sample size:
#' \itemize{
#'   \item \code{"first"}: the smallest \eqn{N}{N} reaching the target power.
#'     The power may fall below the target again for slightly larger \eqn{N}{N}.
#'   \item \code{"stays_above"}: the smallest \eqn{N}{N} from which the power
#'     stays at or above the target for all larger \eqn{N}{N}.
#' }
#'
#' The methodology follows the joint sample size framework for single-arm
#' confirmatory diagnostic accuracy studies described in Stark, M., & Zapf, A. (2020).
#'
#' In the single-test design there is no distinction between superiority and non-inferiority.
#'
#' @references
#' Stark, M., & Zapf, A. (2020).
#' Sample size calculation and re-estimation based on the prevalence
#' in a single-arm confirmatory diagnostic accuracy study.
#' \emph{Statistical Methods in Medical Research}, 29(10), 2958–2971.
#' \doi{10.1177/0962280220913588}
#' @param power Overall desired power
#' @param prev Disease prevalence
#' @param se.e Sensitivity under alternative hypothesis
#' @param se.min Sensitivity under null hypothesis
#' @param sp.e Specificity under alternative hypothesis
#' @param sp.min Specificity under null hypothesis
#' @param alpha Type I error level.
#' @param exact If \code{TRUE}, use the exact binomial power instead of the
#'   normal approximation (see Details). Default \code{FALSE}.
#' @param exact_rule Only for \code{exact = TRUE}: \code{"first"} (default) or
#'   \code{"stays_above"}, see Details.
#' @param ci_method Normal approximation method: "wald" (default) or "logit".
#' @param alternative "one.sided" (default) or "two.sided".
#' @author Ainesh Sewak & Maia Muresan (DCR Bern)
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
#'   \item{target_power}{Requested power (only when \code{power} was given).}
#'   \item{exact_rule}{Rule used (only for \code{exact = TRUE} and a given \code{power}).}
#' }
#' @noRd
#' @examples
#' # Example 1: Sample size calculation for a single diagnostic test
#' # Objective: demonstrate that sensitivity and specificity exceed
#' # pre-specified minimum acceptable values
#' ss_opt_bin(
#'   power = 0.80,        # desired power
#'   prev = 0.5,          # disease prevalence
#'   se.e = 0.85,         # expected sensitivity
#'   se.min = 0.65,       # minimum acceptable sensitivity
#'   sp.e = 0.85,         # expected specificity
#'   sp.min = 0.55,       # minimum acceptable specificity
#'   alpha = 0.025        # one-sided significance level
#' )
#' #same
#' ss_opt_bin(
#'   power = 0.80,        # desired power
#'   prev = 0.5,          # disease prevalence
#'   se.e = 0.85,         # expected sensitivity
#'   se.min = 0.65,       # minimum acceptable sensitivity
#'   sp.e = 0.85,         # expected specificity
#'   sp.min = 0.55,       # minimum acceptable specificity
#'   alpha = 0.025,        # one-sided significance level (default)
#'   ci_method = "wald",   # default
#'   alternative = "one.sided" #default
#' )
#' # Example 2: Same design using logit-based confidence interval
#' ss_opt_bin(
#'   power = 0.80,
#'   prev = 0.5,
#'   se.e = 0.85,
#'   se.min = 0.65,
#'   sp.e = 0.85,
#'   sp.min = 0.55,
#'   alpha = 0.025,
#'   ci_method = "logit",
#'   alternative = "one.sided"
#' )
ss_opt_bin <- function(power,
                       prev,
                       se.e,se.min,
                       sp.e,sp.min,
                       alpha = 0.025, ci_method=c("wald","logit"),
                       alternative=c("one.sided","two.sided"),
                       exact = FALSE,
                       exact_rule = c("first", "stays_above")) {

  ci_method <- match.arg(ci_method)
  alternative<-match.arg(alternative)
  exact_rule <- match.arg(exact_rule)

  # -----------------------------
  # Input validation
  # -----------------------------
  check_inputs(power = power, prev = prev, alpha = alpha)
  check_single(se.e, se.min, sp.e, sp.min)

  # smallest total N such that the joint power >= power
  pwr_at <- function(N, exact = FALSE) {
    power_opt_bin(N = N, prev = prev,
                  se.e = se.e, se.min = se.min,
                  sp.e = sp.e, sp.min = sp.min,
                  alpha = alpha, ci_method = ci_method,
                  alternative = alternative, exact = exact)
  }

  N <- find_n(function(N) pwr_at(N)$power, power)

  if (exact) {
    z_alpha <- z_critical(alpha, alternative)
    exact_pw <- function(N) {
      ex <- exact_power_bin(N, prev, se.e, se.min, sp.e, sp.min,
                            z_alpha, ci_method)
      ex$power_se * ex$power_sp
    }
    N <- exact_n(exact_pw, power, N_asym = N, rule = exact_rule)
  }
  result <- pwr_at(N, exact = exact)

  result$method <- power_method("Single-test design",
                                single_approx(ci_method, exact), "samplesize")
  result$target_power <- power
  if (exact) result$exact_rule <- exact_rule
  result
}

