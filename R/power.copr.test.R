#' Power and Sample Size for Co-Primary Diagnostic Accuracy Endpoints - Single-test design
#'
#' @description
#' Computes either the required total sample size or the achieved overall power
#' of a single-test confirmatory diagnostic accuracy study with co-primary
#' endpoints (sensitivity and specificity), in the style of
#' \code{\link[stats]{power.t.test}}. Exactly one of \code{N} and \code{power}
#' must be \code{NULL}; the one that is \code{NULL} is calculated.
#'
#' Both endpoints are tested for superiority against fixed null values and
#' study success requires rejection of both null hypotheses (intersection-union
#' test). See \code{\link{ss_opt_bin}} and \code{\link{power_opt_bin}} for the
#' methodology.
#'
#' @param N Total sample size, or \code{NULL} to calculate it.
#' @param power Desired overall power, or \code{NULL} to calculate it.
#' @param prev Disease prevalence (0 < prev < 1).
#' @param se.e Sensitivity under the alternative hypothesis.
#' @param se.min Sensitivity under the null hypothesis.
#' @param sp.e Specificity under the alternative hypothesis.
#' @param sp.min Specificity under the null hypothesis.
#' @param alpha Type I error level (default 0.025).
#' @param ci_method Normal approximation method: "wald" (default) or "logit".
#' @param alternative "one.sided" (default) or "two.sided".
#' @param rho Correlation between the sensitivity and specificity test
#'   statistics (default 0), see \code{\link{power_opt_bin}}.
#'
#' @return An object of class \code{"dx_bin"}, see \code{\link{ss_opt_bin}}
#' (when \code{power} is given) and \code{\link{power_opt_bin}} (when \code{N}
#' is given).
#'
#' @seealso \code{\link{ss_opt_bin}}, \code{\link{power_opt_bin}}
#' @author Ainesh Sewak & Maia Muresan (DCR Bern)
#' @export
#' @examples
#' # Sample size for 80% power
#' power.copr.test(power = 0.80, prev = 0.5,
#'                 se.e = 0.85, se.min = 0.65,
#'                 sp.e = 0.85, sp.min = 0.55)
#'
#' # Power for a fixed sample size
#' power.copr.test(N = 122, prev = 0.5,
#'                 se.e = 0.85, se.min = 0.65,
#'                 sp.e = 0.85, sp.min = 0.55,
#'                 ci_method = "logit")
power.copr.test <- function(N = NULL,
                            power = NULL,
                            prev,
                            se.e, se.min,
                            sp.e, sp.min,
                            alpha = 0.025,
                            ci_method = c("wald", "logit"),
                            alternative = c("one.sided", "two.sided"),
                            rho = 0) {

  if (sum(vapply(list(N, power), is.null, logical(1))) != 1) {
    stop("exactly one of 'N' and 'power' must be NULL")
  }

  if (is.null(N)) {
    ss_opt_bin(power = power, prev = prev,
               se.e = se.e, se.min = se.min,
               sp.e = sp.e, sp.min = sp.min,
               alpha = alpha, ci_method = ci_method,
               alternative = alternative, rho = rho)
  } else {
    power_opt_bin(N = N, prev = prev,
                  se.e = se.e, se.min = se.min,
                  sp.e = sp.e, sp.min = sp.min,
                  alpha = alpha, ci_method = ci_method,
                  alternative = alternative, rho = rho)
  }
}
