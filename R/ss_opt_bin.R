#' Optimal Sample Size Calculation for Co-Primary Diagnostic Accuracy Endpoints - Single-test design
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
#' error rate is controlled at level \eqn{\alpha} without multiplicity
#' adjustment.
#'
#' The calculation explicitly accounts for disease prevalence, which
#' determines the expected number of diseased and non-diseased subjects.
#'
#' @details
#' Let \eqn{\pi} denote the assumed disease prevalence. The study population
#' is partitioned into diseased (\eqn{N_{Se} = N\pi}) and non-diseased
#' (\eqn{N_{Sp} = N(1-\pi)}) independent subgroups.
#'
#' When prevalence is very low or very high, one subgroup may require
#' substantially more subjects than the other. A naive calculation (e.g. taking
#' the maximum of two independently calculated sample sizes) may
#' therefore overpower one endpoint while the other determines the
#' required total sample size.
#'
#' This function implements a joint sample size calculation that accounts
#' for both endpoints simultaneously: the required total sample size is the
#' smallest integer \eqn{N} for which the overall power, as computed by
#' \code{\link{power_opt_bin}}, reaches the target. Neither sensitivity nor
#' specificity is treated as the sole driver of the sample size.
#' Power is based on large-sample normal approximations
#' for binomial proportions using either:
#' \itemize{
#'   \item the Wald approximation on the probability scale, or
#'   \item the logit-transformed normal approximation.
#' }
#'
#' Overall power is the bivariate normal probability that both tests reject
#' (\code{mvtnorm::pmvnorm}), with correlation \code{rho} between the two
#' test statistics. With the default \code{rho = 0} this equals
#' \deqn{
#' Power_{overall} = Power_{Se} \times Power_{Sp}.
#' }
#' \eqn{N} is found by bracketing (doubling) and bisection.
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
#' @param ci_method Normal approximation method: "wald" (default) or "logit".
#' @param alternative "one.sided" (default) or "two.sided".
#' @param rho Correlation between the sensitivity and specificity test
#'   statistics (default 0). Sensitivity and specificity are estimated in
#'   disjoint subgroups, so 0 is the natural value; non-zero values are
#'   intended for sensitivity analyses.
#' @author Ainesh Sewak & Maia Muresan (DCR Bern)
#' @importFrom stats approx pnorm qnorm spline
#' @importFrom mvtnorm pmvnorm TVPACK
#' @return An object of class "dx_bin" containing:
#' \itemize{
#'   \item design specifications(design specifications (alpha level, testing framework, CI method)
#'   \item required total sample size,
#'   \item subgroup sample sizes,
#'   \item target overall power and achieved overall power (\code{overall_power})
#'   \item achieved endpoint-specific powers.
#'
#' }
#' @export
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
                       rho = 0) {

  ci_method <- match.arg(ci_method)
  alternative<-match.arg(alternative)

  # -----------------------------
  # Input validation
  # -----------------------------
  if (power <= 0 || power >= 1)stop("power must be between 0 and 1.")
  if (prev <= 0 || prev >= 1)stop("Prevalence must be between 0 and 1.")
  if (se.e <= se.min)stop("se.e must be greater than se.min for superiority testing.")
  if (sp.e <= sp.min)stop("sp.e must be greater than sp.min for superiority testing.")
  if (alpha <= 0 || alpha >= 1)stop("alpha must be between 0 and 1.")

  # smallest total N such that the joint (bivariate normal) power >= power
  pwr_at <- function(N) {
    power_opt_bin(N = N, prev = prev,
                  se.e = se.e, se.min = se.min,
                  sp.e = sp.e, sp.min = sp.min,
                  alpha = alpha, ci_method = ci_method,
                  alternative = alternative, rho = rho)
  }

  N <- find_n(function(N) pwr_at(N)$overall_power, power)
  result <- pwr_at(N)

  result$analysis <- "samplesize"
  result$target_power <- power
  result
}

