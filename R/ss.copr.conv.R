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
#' \right)^2,}{n = ((z_(1-alpha) * sqrt(V0) + z_(1-beta) * sqrt(VA)) / d)^2,}
#' with \eqn{z_{1-\alpha}}{z_(1-alpha)} the one-sided (or two-sided,
#' \eqn{\alpha/2}{alpha/2}) critical value, \eqn{z_{1-\beta}}{z_(1-beta)} the
#' normal quantile of the per-endpoint power \eqn{1-\beta}{1 - beta},
#' \eqn{d}{d} the difference between the assumed value (\code{se.e} or
#' \code{sp.e}) and the minimum (\code{se.min} or \code{sp.min}), and
#' \eqn{V_0}{V0}, \eqn{V_A}{VA} the variances at the minimum and at the
#' assumed value. On the Wald scale \eqn{d}{d} is the plain difference and
#' \eqn{V = p(1-p)}{V = p * (1 - p)}; on the logit scale \eqn{d}{d} is the
#' difference of the logits and \eqn{V = 1/(p(1-p))}{V = 1 / (p * (1 - p))}.
#' The power at a given total sample size does not depend on
#' how that sample size was chosen, so the power of a conventional design is
#' obtained from \code{\link{power.copr.test}}.
#'
#' @inheritParams power.copr.test
#' @param power_endpoint Power used for each endpoint in step 1. The default
#'   \code{1 - (1 - power) / 2} splits the allowed type II error equally over
#'   the two endpoints (0.9 for an overall power of 0.8, as in Stark et al.).
#'   It must lie in \code{[sqrt(power), 1)}: the endpoints are independent, so
#'   if both have power \code{power_endpoint} the overall power is its square,
#'   which reaches the target \code{power} only from \code{sqrt(power)} on.
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
  warn_alpha(alpha, alternative)
  check_single(se.e, se.min, sp.e, sp.min)
  if (power_endpoint < sqrt(power) - 1e-12 || power_endpoint >= 1)
    stop("power_endpoint must be in [sqrt(power), 1): with independent ",
         "endpoints the overall power is the product of the endpoint powers.",
         call. = FALSE)

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
