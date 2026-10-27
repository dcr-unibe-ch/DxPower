#' Power calculation for Single-test
#' Diagnostic Accuracy Study with Co-Primary Endpoints
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
#' at level \eqn{\alpha} without multiplicity adjustment.
#'
#' The calculation explicitly incorporates the assumed disease prevalence,
#' which determines the expected number of diseased and non-diseased subjects.
#'
#' @details
#' Let \eqn{\pi} denote the assumed disease prevalence. The study population
#' is partitioned into diseased (\eqn{N_{Se} = N\pi}) and non-diseased
#' (\eqn{N_{Sp} = N(1-\pi)}) independent subgroups.
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
#' critical values, computed as a bivariate normal probability
#' (\code{mvtnorm::pmvnorm}, TVPACK) with correlation \code{rho}. With the
#' default \code{rho = 0} this equals
#'
#' \deqn{
#' Power_{overall} = Power_{Se} \times Power_{Sp}.
#' }
#'
#' The function returns the achieved endpoint-specific powers and
#' their joint probability for the specified total sample size.
#' @seealso \code{\link{ss_opt_bin}}
#' @references
#' Stark, M., & Zapf, A. (2020). Sample size calculation and re-estimation
#' based on the prevalence in a single-arm confirmatory diagnostic accuracy
#' study. *Statistical Methods in Medical Research*, 29(10), 2958–2971.
#' \doi{10.1177/0962280220913588}
#' The power calculation implemented in this function is consistent with
#' the asymptotic normal framework underlying the sample size methodology
#' described in the above reference.
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
#' @param rho Correlation between the sensitivity and specificity test
#'   statistics (default 0). Sensitivity and specificity are estimated in
#'   disjoint subgroups, so 0 is the natural value; non-zero values are
#'   intended for sensitivity analyses.
#' @author Maia Muresan (DCR Bern)
#'
#' @importFrom stats approx pnorm qnorm spline
#' @importFrom mvtnorm pmvnorm TVPACK
#' @return An object of class `"dx_bin"` containing:
#' \itemize{
#'   \item design specifications (design specifications (alpha level, testing framework, CI method)
#'   \item achieved overall power,
#'   \item endpoint-specific powers,
#'   \item subgroup sample sizes.
#' }
#'
#' @export
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
#'   alpha = 0.025,      # one-sided significance level (default)#'   ci_method = "wald"  # default
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
                          rho = 0
) {

  ci_method   <- match.arg(ci_method)
  alternative <- match.arg(alternative)
  # -----------------------------
  # Input validation
  # -----------------------------
  if (rho <= -1 || rho >= 1) stop("rho must be strictly between -1 and 1.")
  if (N <= 0) stop("N must be positive.")
  if (prev <= 0 || prev >= 1) stop("Prevalence must be between 0 and 1.")
  if (se.e <= se.min) stop("se.e must be greater than se.min.")
  if (sp.e <= sp.min) stop("sp.e must be greater than sp.min.")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1.")

  # -----------------------------
  # Common quantities
  # -----------------------------
  #critical value
  z_alpha <- if (alternative == "one.sided") {
    qnorm(1 - alpha)
  } else {
    qnorm(1 - alpha/2)
  }

  #sample sizes per subgroup
  Nd <- N * prev           # number of diseased
  Nh <- N * (1 - prev)     # number of non-diseased

  logit <- function(p) log(p / (1 - p))


  # -----------------------------
  # Define scale-specific quantities
  # -----------------------------
  if (ci_method == "wald") {

    d_se <- se.e - se.min
    d_sp <- sp.e - sp.min

    v_se.min <- se.min * (1 - se.min)
    v_se.e <- se.e * (1 - se.e)

    v_sp.min <- sp.min * (1 - sp.min)
    v_sp.e <- sp.e * (1 - sp.e)


  } else if (ci_method == "logit") {

    d_se <- logit(se.e) - logit(se.min)
    d_sp <- logit(sp.e) - logit(sp.min)

    v_se.min <- 1 / (se.min * (1 - se.min))
    v_se.e <- 1 / (se.e * (1 - se.e))

    v_sp.min <- 1 / (sp.min * (1 - sp.min))
    v_sp.e <- 1 / (sp.e * (1 - sp.e))
  }

  # -----------------------------
  # Endpoint-specific power
  # -----------------------------
  #Sensitivity power
  z_beta_se <- (sqrt(Nd) * d_se - z_alpha * sqrt(v_se.min)) / sqrt(v_se.e)
  power_se <- pnorm(z_beta_se)

  # Specificity power
  z_beta_sp <- (sqrt(Nh) * d_sp - z_alpha * sqrt(v_sp.min)) / sqrt(v_sp.e)
  power_sp <- pnorm(z_beta_sp)

  # -----------------------------
  # Overall power (IUT framework)
  # -----------------------------
  # P(both tests reject) = P(Z_se < z_beta_se, Z_sp < z_beta_sp) for standard
  # normal Z with correlation rho (rho = 0 gives power_se * power_sp)
  R <- matrix(c(1, rho, rho, 1), 2)
  overall_power <- as.numeric(
    mvtnorm::pmvnorm(lower = c(-Inf, -Inf),
                     upper = c(z_beta_se, z_beta_sp),
                     corr  = R,
                     algorithm = mvtnorm::TVPACK())
  )


  # -----------------------------
  # Return structured object
  # -----------------------------
  result <- list(
    alternative=alternative,
    alpha = alpha,
    hypothesis=list(
      sensitivity = "sup",
      specificity = "sup"
    ),
    design  = "single-test design",
    analysis = "power",
    ci_method = ci_method,
    overall_power = overall_power,
    power_se = power_se,
    power_sp = power_sp,
    rho = rho,
    N_se = Nd,
    N_sp = Nh,
    N_tot = N,
    prev = prev
  )

  class(result) <- "dx_bin"
  return(result)
}


