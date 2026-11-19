#' Empirical (simulated) power - Single-test design
#'
#' @description
#' Monte Carlo estimate of the probability that both co-primary endpoints are
#' declared significant in a single-test design: the observed sensitivity
#' (specificity) is compared with the fixed minimum acceptable value
#' \code{se.min} (\code{sp.min}). The tests are the ones assumed by
#' \code{\link{power.copr.test}}, so the result can be compared directly with
#' the asymptotic power.
#'
#' @details
#' The number of true positives is drawn from
#' \eqn{\mathrm{Bin}(n_{dis}, Se_E)}{Bin(n_dis, se.e)} and the number of true
#' negatives from \eqn{\mathrm{Bin}(n_{nondis}, Sp_E)}{Bin(n_nondis, sp.e)},
#' independently, for all \code{n_sim}
#' replications at once (no loop). For an observed proportion
#' \eqn{\hat p}{p_hat} and null value \eqn{p_0}{p0} an endpoint is rejected if
#' \deqn{\sqrt{n}\,(g(\hat p) - g(p_0)) / \sqrt{V_0} > z_{1-\alpha},}{sqrt(n) * (g(p_hat) - g(p0)) / sqrt(V0) > z_(1-alpha),}
#' with \eqn{g}{g} the identity (\code{"wald"}) or the logit (\code{"logit"}).
#' \eqn{V_0}{V0} is the null variance of \eqn{\sqrt{n}\,g(\hat p)}{sqrt(n) * g(p_hat)}:
#' \eqn{V_0 = p_0(1-p_0)}{V0 = p0 * (1 - p0)} on the Wald scale and
#' \eqn{V_0 = 1/(p_0(1-p_0))}{V0 = 1 / (p0 * (1 - p0))} on the logit scale
#' (the delta-method variance of the logit of a proportion). The Wald version uses the
#' variance under the null and therefore coincides with the Wilson score
#' interval: rejecting is the same as the lower limit of
#' \code{prop.test(x, n, conf.level = 1 - 2 * alpha, correct = FALSE)}
#' exceeding \code{p_0}.
#'
#' @param n_dis,n_nondis Number of diseased and non-diseased individuals.
#' @inheritParams power.copr.test
#' @param n_sim Number of simulation replications.
#'
#' @return Named numeric vector with the empirical \code{overall} power and the
#'   empirical power of the sensitivity (\code{se}) and specificity
#'   (\code{sp}) test. The Monte Carlo standard error of each is at most
#'   \eqn{0.5/\sqrt{n_{sim}}}{0.5 / sqrt(n_sim)}.
#'
#' @seealso \code{\link{power.copr.test}},\code{\link{power.copr.comp.emp}}
#' @importFrom stats rbinom qlogis
#' @export
#' @examples
#' set.seed(1)
#' power.copr.emp(n_dis = 61, n_nondis = 61,
#'                     se.e = 0.85, se.min = 0.65,
#'                     sp.e = 0.85, sp.min = 0.55)
power.copr.emp <- function(n_dis, n_nondis,
                                se.e, se.min,
                                sp.e, sp.min,
                                alpha = 0.025,
                                ci_method = c("wald", "logit"),
                                alternative = c("one.sided", "two.sided"),
                                n_sim = 10000) {

  ci_method   <- match.arg(ci_method)
  alternative <- match.arg(alternative)
  check_num(alpha, "alpha", 0, 1)
  warn_alpha(alpha, alternative)
  check_count(n_dis, "n_dis")
  check_count(n_nondis, "n_nondis")
  check_count(n_sim, "n_sim")
  check_single(se.e, se.min, sp.e, sp.min)

  z <- z_critical(alpha, alternative)

  rej_se <- reject_single(rbinom(n_sim, n_dis,    se.e), n_dis,    se.min, z, ci_method)
  rej_sp <- reject_single(rbinom(n_sim, n_nondis, sp.e), n_nondis, sp.min, z, ci_method)

  c(overall = mean(rej_se & rej_sp), se = mean(rej_se), sp = mean(rej_sp))
}

# Decision rule of the single-test design, given the observed counts `x` out
# of `n` (vectorised over x): TRUE if the Wald (null variance) or logit
# statistic exceeds the critical value `z`. Kept separate from the simulation
# so that it can be tested at known counts.
reject_single <- function(x, n, min, z, ci_method = c("wald", "logit")) {
  ci_method <- match.arg(ci_method)
  p_hat <- x / n
  stat <- if (ci_method == "wald") {
    (p_hat - min) / sqrt(min * (1 - min) / n)
  } else {
    (qlogis(p_hat) - qlogis(min)) * sqrt(n * min * (1 - min))
  }
  stat > z
}

#' Empirical (simulated) power - Unpaired comparative design
#'
#' @description
#' Monte Carlo estimate of the probability that both co-primary endpoints are
#' declared significant in an unpaired comparative design. Unlike the
#' single-test design, the comparator proportion is itself estimated from a
#' sample of the same size, so the result is \emph{not} the same as
#' \code{\link{power.copr.emp}} with \code{se.min = se.c}: the test
#' statistic has the variance of a difference of two estimated proportions.
#'
#' @details
#' For each endpoint, \eqn{\hat p_C \sim \mathrm{Bin}(n, p_C)/n}{p_C_hat ~ Bin(n, p_C) / n} and
#' \eqn{\hat p_E \sim \mathrm{Bin}(n, p_E)/n}{p_E_hat ~ Bin(n, p_E) / n} are drawn independently (\code{n_dis}
#' per arm for the sensitivity, \code{n_nondis} per arm for the specificity)
#' and the Wald statistic
#' \eqn{(\hat p_E - \hat p_C + \delta) / \widehat{\mathrm{se}} > z_{1-\alpha}}{(p_E_hat - p_C_hat + delta) / se_hat > z_(1-alpha)}
#' is evaluated, where \eqn{\delta}{delta} is the non-inferiority margin (0 for
#' superiority) and the standard error \eqn{\widehat{\mathrm{se}}}{se_hat} uses either the unpooled variance or the
#' pooled proportion \eqn{(\hat p_C + \hat p_E)/2}{(p_C_hat + p_E_hat) / 2}. Replications with a zero
#' standard error and a zero numerator count as not rejected.
#'
#' @param n_dis,n_nondis Number of diseased and non-diseased individuals
#'   \emph{per arm}.
#' @inheritParams power.copr.comp.test
#' @param variance Variance estimate in the test statistic: \code{"unpooled"}
#'   (default) or \code{"pooled"}. (The \code{"null"} option of
#'   \code{\link{power.copr.comp.test}} uses the true comparator variance
#'   and has no empirical counterpart.)
#' @param n_sim Number of simulation replications.
#'
#' @return Named numeric vector with the empirical \code{overall}, \code{se}
#'   and \code{sp} power.
#'
#' @seealso \code{\link{power.copr.comp.test}}, \code{\link{power.copr.emp}}
#' @export
#' @examples
#' set.seed(1)
#' power.copr.comp.emp(n_dis = 150, n_nondis = 100,
#'                              se.c = 0.75, se.e = 0.85,
#'                              sp.c = 0.65, sp.e = 0.75,
#'                              alpha = 0.025)
power.copr.comp.emp <- function(n_dis, n_nondis,
                                         se.c, se.e,
                                         sp.c, sp.e,
                                         alpha = 0.025,
                                         delta.se = 0, delta.sp = 0,
                                         type.se = c("sup", "ni"),
                                         type.sp = c("sup", "ni"),
                                         variance = c("unpooled", "pooled"),
                                         alternative = c("one.sided", "two.sided"),
                                         n_sim = 10000) {

  type.se     <- match.arg(type.se)
  type.sp     <- match.arg(type.sp)
  variance    <- match.arg(variance)
  alternative <- match.arg(alternative)
  check_num(alpha, "alpha", 0, 1)
  warn_alpha(alpha, alternative)
  check_count(n_dis, "n_dis")
  check_count(n_nondis, "n_nondis")
  check_count(n_sim, "n_sim")
  for (nm in c("se.c", "se.e", "sp.c", "sp.e")) check_num(get(nm), nm, 0, 1)

  z <- z_critical(alpha, alternative)

  # simulate both arms of one endpoint and apply the decision rule
  reject <- function(n, c, e, delta, type) {
    margin <- if (type == "ni") delta else 0
    reject_unpaired(rbinom(n_sim, n, c), rbinom(n_sim, n, e), n, margin, z,
                    variance)
  }

  rej_se <- reject(n_dis,    se.c, se.e, delta.se, type.se)
  rej_sp <- reject(n_nondis, sp.c, sp.e, delta.sp, type.sp)

  c(overall = mean(rej_se & rej_sp), se = mean(rej_se), sp = mean(rej_sp))
}

# Decision rule of the unpaired design, given the comparator count `x` and the
# experimental count `y` out of `n` per arm (vectorised): TRUE if the Wald
# statistic (y/n - x/n + margin) / se_hat exceeds `z`, with the unpooled or
# pooled variance estimate. A zero standard error with a zero numerator (NaN)
# counts as not rejected. `margin` is delta for "ni" and 0 for "sup".
reject_unpaired <- function(x, y, n, margin, z, variance = c("unpooled", "pooled")) {
  variance <- match.arg(variance)
  x <- x / n
  y <- y / n
  se <- if (variance == "pooled") {
    pb <- (x + y) / 2
    sqrt(2 * pb * (1 - pb) / n)
  } else {
    sqrt((x * (1 - x) + y * (1 - y)) / n)
  }
  rej <- (y - x + margin) / se > z
  rej[is.na(rej)] <- FALSE
  rej
}
