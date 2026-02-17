#' Power calculation for Co-Primary Diagnostic Accuracy Endpoints
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
#' Overall power is therefore computed as:
#'
#' \deqn{
#' Power_{overall} = Power_{Se} \times Power_{Sp}.
#' }
#'
#' The function returns the achieved endpoint-specific powers and
#' their joint probability for the specified total sample size.
#'
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
#' @param se1 Sensitivity under the alternative hypothesis.
#' @param se0 Sensitivity under the null hypothesis.
#' @param sp1 Specificity under the alternative hypothesis.
#' @param sp0 Specificity under the null hypothesis.
#' @param alpha Type I error level.
#' @param ci_method Normal approximation method: "wald" (default) or "logit".
#' @param alternative "one.sided" (default) or "two.sided".
#' @author Maia Muresan (DCR Bern)
#'
#' @importFrom stats approx pnorm qnorm spline
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
#' power_opt_bin(N=122,prev=0.5,se1=0.85, se0=0.65,sp1=0.85, sp0=0.55,alpha = 0.025)
#' power_opt_bin(N=122,prev=0.5,se1=0.85, se0=0.65,sp1=0.85, sp0=0.55,alpha = 0.025, ci_method="logit", alternative="two.sided")

power_opt_bin <- function(N,
                          prev,
                          se1, se0,
                          sp1, sp0,
                          alpha = 0.025, ci_method=c("wald", "logit"),
                          alternative=c("one.sided","two.sided")
) {

  ci_method   <- match.arg(ci_method)
  alternative <- match.arg(alternative)
  # -----------------------------
  # Input validation
  # -----------------------------
  if (N <= 0) stop("N must be positive.")
  if (prev <= 0 || prev >= 1) stop("Prevalence must be between 0 and 1.")
  if (se1 <= se0) stop("se1 must be greater than se0.")
  if (sp1 <= sp0) stop("sp1 must be greater than sp0.")
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

    delta_se <- se1 - se0
    delta_sp <- sp1 - sp0

    se_se0 <- sqrt(se0 * (1 - se0))
    se_se1 <- sqrt(se1 * (1 - se1))

    se_sp0 <- sqrt(sp0 * (1 - sp0))
    se_sp1 <- sqrt(sp1 * (1 - sp1))


  } else if (ci_method == "logit") {

    delta_se <- logit(se1) - logit(se0)
    delta_sp <- logit(sp1) - logit(sp0)

    se_se0 <- sqrt(1 / (se0 * (1 - se0)))
    se_se1 <- sqrt(1 / (se1 * (1 - se1)))

    se_sp0 <- sqrt(1 / (sp0 * (1 - sp0)))
    se_sp1 <- sqrt(1 / (sp1 * (1 - sp1)))
  }

  # -----------------------------
  # Endpoint-specific power
  # -----------------------------
  #Sensitivity power
  z_beta_se <- (sqrt(Nd) * delta_se - z_alpha * se_se0) / se_se1
  pwr_se <- pnorm(z_beta_se)

  # Specificity power
  z_beta_sp <- (sqrt(Nh) * delta_sp - z_alpha * se_sp0) / se_sp1
  pwr_sp <- pnorm(z_beta_sp)

  # -----------------------------
  # Overall power (IUT framework)
  # -----------------------------

  overall <- pwr_se * pwr_sp


  # -----------------------------
  # Return structured object
  # -----------------------------
  result <- list(
    alternative=alternative,
    alpha = alpha,
    hypothesis="superiority",
    framework = "co-primary (IUT)",
    type = "power",
    ci_method = ci_method,
    overall = overall,
    pwr_se = pwr_se,
    pwr_sp = pwr_sp,
    N_se = Nd,
    N_sp = Nh,
    N_tot = N,
    prev = prev
  )

  class(result) <- "dx_bin"
  return(result)
}

