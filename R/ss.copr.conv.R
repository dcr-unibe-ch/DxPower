#' Conventional sample size calculation for co-primary diagnostic accuracy endpoints - Single-test design
#'
#' @description
#' Sample size by the conventional three-step approach (Stark et al., 2022):
#' (1) calculate the number of diseased individuals needed for the
#' sensitivity and the number of non-diseased individuals needed for the
#' specificity, each with its own endpoint power, (2) divide these numbers by
#' the prevalence and by one minus the prevalence, respectively, to obtain the
#' total sample size needed for each endpoint, and (3) take the maximum.
#'
#' The endpoint with the smaller requirement is then overpowered, so the
#' overall power is typically above the target, in particular for unbalanced
#' prevalence. Compare with \code{\link{power.copr.test}}, which finds the smallest
#' total sample size reaching the target overall power.
#'
#' @details
#' The per-endpoint sample sizes use the same large-sample approximations as
#' \code{\link{power.copr.test}} (Wald or logit scale), so
#' \deqn{n = \left(\frac{z_{1-\alpha}\sqrt{V_0} + z_{1-\beta}\sqrt{V_A}}{d}
#' \right)^2,}
#' with \eqn{z_{1-\alpha}} the one-sided (or two-sided, \eqn{\alpha/2})
#' critical value. The power at a given total sample size does not depend on
#' how that sample size was chosen, so the power of a conventional design is
#' obtained from \code{\link{power.copr.test}}.
#'
#' @inheritParams power.copr.test
#' @param power_endpoint Power used for each endpoint in step 1. The default
#'   \code{1 - (1 - power) / 2} splits the allowed type II error equally over
#'   the two endpoints (0.9 for an overall power of 0.8, as in Stark et al.).
#'
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
#'   \item{N_se_required}{Total sample size needed for the sensitivity endpoint.}
#'   \item{N_sp_required}{Total sample size needed for the specificity endpoint.}
#'   \item{power_endpoint}{Power used for each endpoint in step 1.}
#' }
#'
#' @seealso \code{\link{power.copr.test}}
#' @references
#' Stark, M., Hesse, M., Brannath, W., & Zapf, A. (2022). Blinded sample size
#' re-estimation in a comparative diagnostic accuracy study.
#' \emph{BMC Medical Research Methodology}, 22, 115.
#' \doi{10.1186/s12874-022-01564-2}
#' @importFrom stats qlogis
#' @export
#' @examples
#' ss.copr.conv(power = 0.80, prev = 0.5,
#'             se.e = 0.85, se.min = 0.65,
#'             sp.e = 0.85, sp.min = 0.55)
#'
#' # unbalanced prevalence: the specificity is overpowered
#' ss.copr.conv(power = 0.80, prev = 0.2,
#'             se.e = 0.85, se.min = 0.65,
#'             sp.e = 0.85, sp.min = 0.55)
ss.copr.conv <- function(power,
                        prev,
                        se.e, se.min,
                        sp.e, sp.min,
                        alpha = 0.025, ci_method = c("wald", "logit"),
                        alternative = c("one.sided", "two.sided"),
                        power_endpoint = 1 - (1 - power) / 2) {

  ci_method   <- match.arg(ci_method)
  alternative <- match.arg(alternative)

  check_inputs(power = power, prev = prev, alpha = alpha)
  check_single(se.e, se.min, sp.e, sp.min)
  if (power_endpoint < power || power_endpoint >= 1)
    stop("power_endpoint must be in [power, 1).")

  z_alpha <- z_critical(alpha, alternative)
  z_beta  <- qnorm(power_endpoint)

  # step 1: diseased / non-diseased needed for each endpoint
  n_needed <- function(e, min) {
    s <- endpoint_scale(e, min, ci_method)
    ((z_alpha * sqrt(s$v_min) + z_beta * sqrt(s$v_e)) / s$d)^2
  }
  n_se <- ceiling(n_needed(se.e, se.min))
  n_sp <- ceiling(n_needed(sp.e, sp.min))

  # step 2: total sample size needed for each endpoint
  N_se_required <- ceiling(n_se / prev)
  N_sp_required <- ceiling(n_sp / (1 - prev))

  # step 3: maximum
  N <- max(N_se_required, N_sp_required)

  result <- power_opt_bin(N = N, prev = prev,
                          se.e = se.e, se.min = se.min,
                          sp.e = sp.e, sp.min = sp.min,
                          alpha = alpha, ci_method = ci_method,
                          alternative = alternative)

  result$method         <- power_method(
    "Single-test design",
    paste0(if (ci_method == "wald") "Wald" else "Logit",
           " normal approximation (conventional per-endpoint sample size)"),
    "samplesize")
  result$target_power   <- power
  result$power_endpoint <- power_endpoint
  result$N_se_required  <- N_se_required
  result$N_sp_required  <- N_sp_required
  result
}
