#' Author: Maia Muresan (DCR Bern)


#' Overall Power for Confirmatory Diagnostic Accuracy Study
#'
#' Computes the overall power for co-primary superiority tests on
#' sensitivity and specificity using a normal approximation.
#'
#' The function assumes one-sided hypothesis tests for:
#'   H0: Se <= Se0  vs  H1: Se > Se0
#'   H0: Sp <= Sp0  vs  H1: Sp > Sp0
#'
#' Overall power is calculated under an intersection-union framework,
#' i.e., both tests must reject the null hypothesis.
#'
#' @param N Total sample size
#' @param prev Disease prevalence (0 < prev < 1)
#' @param se1 Sensitivity under the alternative hypothesis
#' @param se0 Sensitivity under the null hypothesis
#' @param sp1 Specificity under the alternative hypothesis
#' @param sp0 Specificity under the null hypothesis
#' @param alpha One-sided type I error (default 0.025)
#'
#' @return A named list containing:
#'   \item{overall}{Overall power}
#'   \item{pwr_se}{Power for sensitivity}
#'   \item{pwr_sp}{Power for specificity}
#'
#' @export



power_opt_bin <- function(N,
                       prev,
                       se1, se0,
                       sp1, sp0,
                       alpha = 0.025) {

  # -----------------------------
  # Input validation
  # -----------------------------

  if (N <= 0)
    stop("N must be positive.")

  if (prev <= 0 || prev >= 1)
    stop("Prevalence must be between 0 and 1.")

  if (se1 <= se0)
    stop("se1 must be greater than se0 for superiority testing.")

  if (sp1 <= sp0)
    stop("sp1 must be greater than sp0 for superiority testing.")

  if (alpha <= 0 || alpha >= 1)
    stop("alpha must be between 0 and 1.")

  # -----------------------------
  # Common quantities
  # -----------------------------

  z_alpha <- qnorm(1 - alpha)

  Nd <- N * prev           # number of diseased
  Nh <- N * (1 - prev)     # number of non-diseased

  # -----------------------------
  # Sensitivity power
  # -----------------------------

  sd_se0 <- sqrt(se0 * (1 - se0))
  sd_se1 <- sqrt(se1 * (1 - se1))
  delta_se <- se1 - se0

  z_beta_se <- (sqrt(Nd) * delta_se - z_alpha * sd_se0) / sd_se1
  pwr_se <- pnorm(z_beta_se)

  # -----------------------------
  # Specificity power
  # -----------------------------

  sd_sp0 <- sqrt(sp0 * (1 - sp0))
  sd_sp1 <- sqrt(sp1 * (1 - sp1))
  delta_sp <- sp1 - sp0

  z_beta_sp <- (sqrt(Nh) * delta_sp - z_alpha * sd_sp0) / sd_sp1
  pwr_sp <- pnorm(z_beta_sp)

  # -----------------------------
  # Overall power (IUT framework)
  # -----------------------------

  overall <- pwr_se * pwr_sp

  return(list(
    overall = overall,
    pwr_se = pwr_se,
    pwr_sp = pwr_sp
  ))
}




