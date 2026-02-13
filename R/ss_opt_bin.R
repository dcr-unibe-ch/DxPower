#' Author: Aiensh Sewak (DCR Bern)

#' Optimal Power Allocation for Sensitivity and Specificity
#'
#' Computes optimal allocation of power between sensitivity and specificity
#' tests such that required sample sizes are equal.
#'
#' @param se1 Sensitivity under alternative hypothesis
#' @param se0 Sensitivity under null hypothesis
#' @param sp1 Specificity under alternative hypothesis
#' @param sp0 Specificity under null hypothesis
#' @param alpha One-sided type I error
#' @param prev Disease prevalence
#' @param pwr Overall desired power
#'
#' @return A list containing allocated powers and required sample sizes
#' @export



ss_opt_bin <- function(se1, se0, sp1, sp0, alpha = 0.025, prev = 0.5, pwr = 0.8) {
  # Fixed values
  qz <- qnorm(alpha, lower.tail = FALSE)
  sd_se0 <- sqrt(se0 * (1 - se0))
  sd_se1 <- sqrt(se1 * (1 - se1))
  sd_sp0 <- sqrt(sp0 * (1 - sp0))
  sd_sp1 <- sqrt(sp1 * (1 - sp1))

  # Define function for root finding
  root_fun <- function(pwr_se) {
    beta_sp <- (pwr_se - pwr) / pwr_se
    pwr_sp  <- 1 - beta_sp
    qpwr_se <- qnorm(pwr_se)
    qpwr_sp <- qnorm(pwr_sp)

    n_se <- (qz * sd_se0 + qpwr_se * sd_se1)^2 / (prev * (se0 - se1)^2)
    n_sp <- (qz * sd_sp0 + qpwr_sp * sd_sp1)^2 / ((1 - prev) * (sp0 - sp1)^2)
    return(n_se - n_sp)
  }

  # Search for root using splines
  eps <- 1e-6
  grd <- seq(pwr + eps , 1 - eps, length.out = 100000)
  s <- spline(grd, root_fun(grd))
  pwr_se <- approx(x = s$y, y = s$x, xout = 0)$y
  beta_sp <- (pwr_se - pwr) / pwr_se
  pwr_sp <- 1 - beta_sp

  # Compute z-scores
  qpwr_se <- qnorm(pwr_se)
  qpwr_sp <- qnorm(pwr_sp)

  # Sample sizes without prevalence
  N_se <- (qz * sd_se0 + qpwr_se * sd_se1)^2 / (se0 - se1)^2
  N_sp <- (qz * sd_sp0 + qpwr_sp * sd_sp1)^2 / (sp0 - sp1)^2

  return(list(
    pwr_se = pwr_se,
    pwr_sp = pwr_sp,
    N_se = N_se,
    N_sp = N_sp,
    N_tot = N_se + N_sp
  ))
}
