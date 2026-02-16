#' Optimal Sample Size Determination for Co-Primary Diagnostic Accuracy Endpoints
#'
#' Computes the required total sample size for a confirmatory diagnostic
#' accuracy study with co-primary endpoints (sensitivity and specificity)
#' under a superiority hypothesis framework.
#'
#' The function determines the optimal allocation of endpoint-specific
#' powers such that both sensitivity and specificity require the same
#' total sample size while achieving a prespecified overall power.
#'
#' Overall power is defined under an intersection–union test (IUT) framework,
#' i.e., both null hypotheses must be rejected:
#'
#'#' Sample size formulas are based on large-sample normal approximations
#' using either the Wald scale or the logit-transformed scale.
#'
#' @param pwr Overall desired power
#' @param prev Disease prevalence
#' @param se1 Sensitivity under alternative hypothesis
#' @param se0 Sensitivity under null hypothesis
#' @param sp1 Specificity under alternative hypothesis
#' @param sp0 Specificity under null hypothesis
#' @param alpha Type I error level.
#' @param ci_method Normal approximation method: "wald" (default) or "logit".
#' @param alternative "one.sided" (default) or "two.sided".
#' @author Ainesh Sewak & Maia Muresan(DCR Bern)
#' @return An object of class "dx_bin".
#' @export

ss_opt_bin <- function(pwr,
                       prev,
                       se1,se0,
                       sp1,sp0,
                       alpha = 0.025, ci_method=c("wald","logit"),
                       alternative=c("one.sided","two.sided")) {

  ci_method <- match.arg(ci_method)
  alternative<-match.arg(alternative)

  # -----------------------------
  # Input validation
  # -----------------------------
  if (pwr <= 0 || pwr >= 1)stop("pwr must be between 0 and 1.")
  if (prev <= 0 || prev >= 1)stop("Prevalence must be between 0 and 1.")
  if (se1 <= se0)stop("se1 must be greater than se0 for superiority testing.")
  if (sp1 <= sp0)stop("sp1 must be greater than sp0 for superiority testing.")
  if (alpha <= 0 || alpha >= 1)stop("alpha must be between 0 and 1.")

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

  logit <- function(p) log(p / (1 - p))
  # -----------------------------
  # Define scale-specific quantities
  # -----------------------------
  if (ci_method == "wald") {

    delta_se <- se1 - se0
    delta_sp <- sp1 - sp0

    v_se0 <- se0 * (1 - se0)
    v_se1 <- se1 * (1 - se1)

    v_sp0 <- sp0 * (1 - sp0)
    v_sp1 <- sp1 * (1 - sp1)

  } else if (ci_method == "logit") {

    delta_se <- logit(se1) - logit(se0)
    delta_sp <- logit(sp1) - logit(sp0)

    v_se0 <- 1 / (se0 * (1 - se0))
    v_se1 <- 1 / (se1 * (1 - se1))

    v_sp0 <- 1 / (sp0 * (1 - sp0))
    v_sp1 <- 1 / (sp1 * (1 - sp1))

  }

  # Define function for root finding
  root_fun <- function(pwr_se) {
    beta_sp <- (pwr_se - pwr) / pwr_se
    pwr_sp  <- 1 - beta_sp
    qpwr_se <- qnorm(pwr_se)
    qpwr_sp <- qnorm(pwr_sp)
    n_se <- (z_alpha * sqrt(v_se0) +
               qpwr_se * sqrt(v_se1))^2 /
      (prev * delta_se^2)

    n_sp <- (z_alpha * sqrt(v_sp0) +
               qpwr_sp * sqrt(v_sp1))^2 /
      ((1 - prev) * delta_sp^2)
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
  N_se <- (z_alpha * sqrt(v_se0) +
             qpwr_se * sqrt(v_se1))^2 /
    (delta_se^2)

  N_sp <- (z_alpha * sqrt(v_sp0) +
             qpwr_sp * sqrt(v_sp1))^2 /
    (delta_sp^2)

  result <- list(
    alternative=alternative,
    alpha = alpha,
    framework = "co-primary (IUT)",
    type = "samplesize",
    ci_method = ci_method,
    target_power = pwr,
    pwr_se = pwr_se,
    pwr_sp = pwr_sp,
    N_se = N_se,
    N_sp = N_sp,
    N_tot = N_se + N_sp,
    prev = prev
  )

  class(result) <- "dx_bin"
  return(result)
}
