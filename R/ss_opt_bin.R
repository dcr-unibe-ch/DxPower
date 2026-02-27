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
#' This function implements an optimal joint sample size calculation
#' that balances both endpoints simultaneously. Endpoint-specific powers
#' are allocated such that:
#'
#' \itemize{
#'   \item the overall power equals the target value,
#'   \item both endpoints imply the same required total sample size,
#'   \item neither sensitivity nor specificity is unnecessarily overpowered.
#' }
#' Sample size formulas are based on large-sample normal approximations
#' for binomial proportions using either:
#' \itemize{
#'   \item the Wald approximation on the probability scale, or
#'   \item the logit-transformed normal approximation.
#' }
#'
#' The endpoint-specific powers are determined numerically by solving
#' a root-finding problem under the overall power constraint:
#' \deqn{
#' Power_{overall} = Power_{Se} \times Power_{Sp}.
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
#' @param ci_method Normal approximation method: "wald" (default) or "logit".
#' @param alternative "one.sided" (default) or "two.sided".
#' @author Ainesh Sewak & Maia Muresan (DCR Bern)
#' @importFrom stats approx pnorm qnorm spline
#' @return An object of class "dx_bin" containing:
#' \itemize{
#'   \item design specifications(design specifications (alpha level, testing framework, CI method)
#'   \item required total sample size,
#'   \item subgroup sample sizes,
#'   \item target overall power
#'   \item allocated endpoint-specific powers.
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
                       alternative=c("one.sided","two.sided")) {

  ci_method <- match.arg(ci_method)
  alternative<-match.arg(alternative)

  # -----------------------------
  # Input validation
  # -----------------------------
  if (power <= 0 || power >= 1)stop("pwr must be between 0 and 1.")
  if (prev <= 0 || prev >= 1)stop("Prevalence must be between 0 and 1.")
  if (se.e <= se.min)stop("se.e must be greater than se.min for superiority testing.")
  if (sp.e <= sp.min)stop("sp.e must be greater than sp.min for superiority testing.")
  if (alpha <= 0 || alpha >= 1)stop("alpha must be between 0 and 1.")

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

  # Define function for root finding
  root_fun <- function(power_se) {
    beta_sp <- (power_se - power) / power_se
    power_sp  <- 1 - beta_sp
    qpower_se <- qnorm(power_se)
    qpower_sp <- qnorm(power_sp)
    n_se <- (z_alpha * sqrt(v_se.min) +
               qpower_se * sqrt(v_se.e))^2 /
      (prev * d_se^2)

    n_sp <- (z_alpha * sqrt(v_sp.min) +
               qpower_sp * sqrt(v_sp.e))^2 /
      ((1 - prev) * d_sp^2)
    return(n_se - n_sp)
  }

  # Search for root using splines
  eps <- 1e-6
  grd <- seq(power + eps , 1 - eps, length.out = 100000)
  s <- spline(grd, root_fun(grd))
  power_se <- approx(x = s$y, y = s$x, xout = 0)$y
  beta_sp <- (power_se - power) / power_se
  power_sp <- 1 - beta_sp

  # Compute z-scores
  qpower_se <- qnorm(power_se)
  qpower_sp <- qnorm(power_sp)

  # Sample sizes without prevalence
  N_se <- (z_alpha * sqrt(v_se.min) +
             qpower_se * sqrt(v_se.e))^2 /
    (d_se^2)

  N_sp <- (z_alpha * sqrt(v_sp.min) +
             qpower_sp * sqrt(v_sp.e))^2 /
    (d_sp^2)

  result <- list(
    alternative=alternative,
    alpha = alpha,
    hypothesis=list(
      sensitivity = "sup",
      specificity = "sup"
    ),
    design  = "single-test design",
    analysis = "samplesize",
    ci_method = ci_method,
    target_power = power,
    power_se = power_se,
    power_sp = power_sp,
    N_se = N_se,
    N_sp = N_sp,
    N_tot = N_se + N_sp,
    prev = prev
  )

  class(result) <- "dx_bin"
  return(result)
}

