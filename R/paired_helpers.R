#' Convert between correlation and discordance for paired binary data
#'
#' @description
#' For one endpoint of a paired design, the two tests give a pair of binary
#' results with marginal probabilities \code{p1} (comparator) and \code{p2}
#' (experimental). The joint distribution is fixed by the marginals and either
#' the discordance \eqn{\psi = P(\mathrm{results\ differ})}{psi = P(results differ)}
#' or the correlation (phi coefficient) \eqn{\rho}{rho}, where
#' \eqn{p_{11}}{p11} is the probability that both tests are positive:
#' \deqn{p_{11} = \frac{p_1 + p_2 - \psi}{2}, \qquad
#'   \rho = \frac{p_{11} - p_1 p_2}{\sqrt{p_1(1-p_1)p_2(1-p_2)}}.}{p11 = (p1 + p2 - psi) / 2,  rho = (p11 - p1 * p2) / sqrt(p1 * (1 - p1) * p2 * (1 - p2)).}
#' \code{rho_to_psi()} and \code{psi_to_rho()} convert between the two.
#' Both parameters are restricted to the range of Stark et al. (2022, Eq. 10),
#' \eqn{|p_1 - p_2| \le \psi \le p_1 + p_2 - 2 p_1 p_2}{|p1 - p2| <= psi <= p1 + p2 - 2 * p1 * p2}.
#' The upper end is the
#' discordance of uncorrelated tests, so this range is
#' \eqn{0 \le \rho \le \rho_{max}}{0 <= rho <= rho_max}: negative correlation
#' (the tests tend to err
#' on different subjects) is not considered. A value outside the range is an
#' error.
#'
#' @param p1,p2 Marginal probabilities (sensitivity or specificity) of the two
#'   tests, strictly between 0 and 1.
#' @param rho Correlation between the two binary results.
#' @param psi Discordance probability.
#' @param name Name used in error messages.
#' @return A single number: \eqn{\psi}{psi} for \code{rho_to_psi()},
#'   \eqn{\rho}{rho} for \code{psi_to_rho()}.
#' @seealso \code{\link{paired_bounds}} for the feasible ranges of an entire
#'   design.
#' @export
#' @examples
#' rho_to_psi(0.81, 0.90, rho = 0.4)
#' psi_to_rho(0.81, 0.90, psi = 0.09)
rho_to_psi <- function(p1, p2, rho, name = "rho") {
  check_num(p1, "p1", 0, 1)
  check_num(p2, "p2", 0, 1)
  check_scalar(rho, name)
  r <- rho_range(p1, p2)
  if (rho < r[1] - 1e-10 || rho > r[2] + 1e-10) {
    stop(name, " = ", rho, " is outside the range of Stark et al. (2022) implied ",
         "by the marginals ", p1, " and ", p2, ": ", round(r[1], 4), " <= ", name,
         " <= ", round(r[2], 4), ".", call. = FALSE)
  }
  p11 <- p1 * p2 + rho * sqrt(p1 * (1 - p1) * p2 * (1 - p2))
  psi <- p1 + p2 - 2 * p11
  pr <- psi_range(p1, p2)
  min(max(psi, pr[1]), pr[2])  # remove floating-point noise at the bounds
}

#' @rdname rho_to_psi
#' @export
psi_to_rho <- function(p1, p2, psi, name = "psi") {
  check_num(p1, "p1", 0, 1)
  check_num(p2, "p2", 0, 1)
  check_scalar(psi, name)
  pr <- psi_range(p1, p2)
  if (psi < pr[1] - 1e-10 || psi > pr[2] + 1e-10) {
    stop(name, " = ", psi, " is outside the range of Stark et al. (2022) implied ",
         "by the marginals ", p1, " and ", p2, ": ", round(pr[1], 4), " <= ", name,
         " <= ", round(pr[2], 4), ".", call. = FALSE)
  }
  psi_to_rho_raw(p1, p2, psi)
}

#' Feasible ranges of discordance and correlation for a paired design
#'
#' @description
#' Before running \code{\link{power.copr.comp.test}} with
#' \code{design = "paired"}, shows which discordance (\code{psi.d},
#' \code{psi.nd}) or correlation (\code{rho.d}, \code{rho.nd}) values are
#' possible given the assumed sensitivities and specificities. The range is
#' that of Stark et al. (2022, Eq. 10),
#' \eqn{|c - e| \le \psi \le c + e - 2ce}{|c - e| <= psi <= c + e - 2 * c * e}
#' (\eqn{c}{c}, \eqn{e}{e}: comparator and experimental sensitivity or
#' specificity),
#' which corresponds to correlations between 0 and a maximum, see
#' \code{\link{rho_to_psi}}. Values outside it are an error in
#' \code{power.copr.comp.test()}. If no discordance or correlation is given
#' there, the lower bound \code{psi.min} is used.
#'
#' @param se.c,se.e Sensitivity of the comparator and the experimental test.
#' @param sp.c,sp.e Specificity of the comparator and the experimental test.
#' @return A data frame with one row per endpoint (the diseased group for
#'   sensitivity, the non-diseased group for specificity) and columns
#'   \code{psi.min}, \code{psi.max} (range of the discordance) and
#'   \code{rho.min}, \code{rho.max} (the correlation range it corresponds to).
#'   \code{rho.min} is always 0, and \code{psi.min} corresponds to
#'   \code{rho.max}.
#' @references
#' Stark, M., Hesse, M., Brannath, W., & Zapf, A. (2022). Blinded sample size
#' re-estimation in a comparative diagnostic accuracy study.
#' \emph{BMC Medical Research Methodology}, 22, 115.
#' \doi{10.1186/s12874-022-01564-2}
#' @seealso \code{\link{rho_to_psi}}, \code{\link{power.copr.comp.test}}
#' @export
#' @examples
#' paired_bounds(se.c = 0.81, se.e = 0.90, sp.c = 0.66, sp.e = 0.80)
paired_bounds <- function(se.c, se.e, sp.c, sp.e) {
  for (nm in c("se.c", "se.e", "sp.c", "sp.e")) check_num(get(nm), nm, 0, 1)
  one <- function(p1, p2) {
    pr <- psi_range(p1, p2)
    rr <- rho_range(p1, p2)
    c(psi.min = pr[1], psi.max = pr[2], rho.min = rr[1], rho.max = rr[2])
  }
  data.frame(endpoint = c("sensitivity", "specificity"),
             rbind(one(se.c, se.e), one(sp.c, sp.e)),
             row.names = NULL)
}
