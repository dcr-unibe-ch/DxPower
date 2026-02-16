#' Overall Power for Co-Primary Diagnostic Accuracy Endpoints
#'
#' Computes the achieved overall power for a confirmatory diagnostic
#' accuracy study with co-primary endpoints (sensitivity and specificity)
#' under a superiority hypothesis framework.
#'
#' Both null hypotheses must be rejected according to an
#' intersection–union test (IUT) principle. The overall power is defined as
#' the joint probability of rejecting both null hypotheses.
#'
#' Endpoint-specific powers are computed using large-sample normal
#' approximations on either the Wald scale or the logit scale.
#'
#'
#' @details
#' The study population is partitioned into diseased and non-diseased
#' subgroups according to the assumed prevalence. Under binomial sampling,
#' the estimators of sensitivity and specificity are asymptotically
#' independent, which justifies the product formulation of overall power.
#'
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
#' @return An object of class "dx_bin".
#' @export


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

  if (alternative == "two.sided" && alpha < 0.05) {
    warning(
      "Two-sided alpha < 0.05 is conservative. ",
      "Typical confirmatory superiority trials use two-sided alpha = 0.05 ",
      "(equivalent to one-sided 0.025)."
    )
  }

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
