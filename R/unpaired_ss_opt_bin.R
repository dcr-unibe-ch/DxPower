#' Optimal Sample Size Calculation for an Unpaired Comparative
#' Diagnostic Accuracy Study with Co-Primary Endpoints
#'
#' @description
#' The function determines the minimal total sample size per group for an unpaired
#' comparative diagnostic accuracy study with co-primary endpoints
#' (sensitivity and specificity) using the intersection–union testing
#' framework (IUT).
#'
#' In the unpaired comparative design, subjects are randomly assigned
#' to either the experimental or the comparator test.All subjects
#' receive by the reference standard.
#'
#' The study is evaluated under the intersection–union testing (IUT)
#' framework: both endpoint-specific null hypotheses must be rejected
#' for overall study success. Under this framework, the global Type I
#' error rate is controlled at level \eqn{\alpha} without multiplicity
#' adjustment.
#'
#' Superiority and non-inferiority hypotheses may be specified
#' independently for sensitivity and specificity.
#'
#' @details
#'
#' Let \eqn{\pi} denote the assumed disease prevalence.
#' With total sample size per group \eqn{N}, the expected subgroup sizes are:
#'
#' \deqn{
#' N_{Se} = N\pi, \qquad N_{Sp} = N(1-\pi).
#' }
#'
#' Sensitivity is evaluated in the diseased subgroup and specificity
#' in the non-diseased subgroup. Since these populations are independent,
#' the overall power under the IUT framework equals:
#'
#' \deqn{
#' Power_{overall} = Power_{Se} \times Power_{Sp}.
#' }
#'
#'
#' Let \eqn{\theta_{Se,E}} and \eqn{\theta_{Se,C}} denote the sensitivities
#' of the experimental and comparator tests, respectively, and
#' \eqn{\theta_{Sp,E}} and \eqn{\theta_{Sp,C}} the corresponding specificities.
#'
#' The global null hypothesis is defined as the union of the
#' endpoint-specific null hypotheses.
#'
#' \strong{Case 1: Superiority for both endpoints}
#'
#' \deqn{
#' H_{0,global} =
#' \left( \theta_{Se,E} \le \theta_{Se,C} \right)
#' \cup
#' \left( \theta_{Sp,E} \le \theta_{Sp,C} \right).
#' }
#'
#' \strong{Case 2: Sensitivity superiority and specificity non-inferiority}
#'
#' Let \eqn{\delta_{Sp} > 0} denote the non-inferiority margin for specificity.
#'
#' \deqn{
#' H_{0,global} =
#' \left( \theta_{Se,E} \le \theta_{Se,C} \right)
#' \cup
#' \left( \theta_{Sp,E} \le \theta_{Sp,C} - \delta_{Sp} \right).
#' }
#'
#' Similar formulations apply for:specificity superiority with sensitivity non-inferiority,
#' or both endpoints non-inferiority.
#'
#' Under IUT, the global null hypothesis is rejected only if both
#' endpoint-specific null hypotheses are rejected.
#'
#' \strong{Sample size formulas}
#'
#' Sample size formulas are based on large-sample normal approximations
#' for the difference of two independent binomial proportions.
#'
#' For superiority hypotheses, the Wald-type approximation is used:
#'
#' \deqn{
#' n =
#' \frac{
#' \left(z_{\alpha^*}\sqrt{V_0} + z_{\beta}\sqrt{V_A}\right)^2
#' }{
#' \Delta^2
#' }
#' }
#' For non-inferiority hypotheses:
#'
#' \deqn{
#' n =
#' \frac{
#' (z_{\alpha} + z_{\beta})^2 V_A
#' }{
#' (\Delta - \delta)^2
#' }
#' }
#'
#' where \eqn{\Delta = \theta_{E} - \theta_{C}} denotes the assumed difference in diagnostic accuracy between comparator and experimental test.
#' and \eqn{\delta} denotes the non-inferiority margin (applicable only for non-inferiority test, otherwise for superiority test the default is 0).
#'
#' The \eqn{V_0} and \eqn{V_A} represent the variances of the
#' difference between \eqn{\theta_{C}} and \eqn{\theta_{E}}
#' under the null (\eqn{V_0}) and alternative hypotheses (\eqn{V_A}), respectively.
#'
#' In the unpaired design, sensitivity and specificity are estimated
#' in independent groups. The variance of the difference between two
#' independent proportions is given by:
#'
#' \deqn{
#' V(\theta_{C} - \theta_{E})
#' =
#' \theta_{C}(1 - \theta_{C})
#' +
#' \theta_{E}(1 - \theta_{E}).
#' }
#'
#' Endpoint-specific powers are allocated numerically such that:
#' \itemize{
#'   \item the overall power equals the prespecified target value,
#'   \item both endpoints imply the same required total sample size
#'   \item neither sensitivity nor specificity is unnecessarily overpowered.
#' }
#' @references
#' Stark, M., Hesse, M., Brannath, W., & Zapf, A. (2022).
#' Blinded sample size re-estimation in a comparative diagnostic
#' accuracy study.
#' \emph{BMC Medical Research Methodology}, 22, 115.
#' \doi{10.1186/s12874-022-01564-2}
#' @param power Overall desired power for the co-primary hypothesis.
#' @param prev Disease prevalence in the study population.
#' @param se.c Sensitivity of the comparator test.
#' @param se.e Sensitivity of the experimental test.
#' @param sp.c Specificity of the comparator test.
#' @param sp.e Specificity of the experimental test.
#' @param alpha Type I error rate per endpoint.
#' @param delta.se Non-inferiority margin for sensitivity
#'   (used only if \code{type.se = "ni"}).
#' @param delta.sp Non-inferiority margin for specificity
#'   (used only if \code{type.sp = "ni"}).
#' @param type.se Type of hypothesis for sensitivity:
#'   \code{"sup"} (superiority, default) or \code{"ni"} (non-inferiority).
#' @param type.sp Type of hypothesis for specificity:
#'   \code{"sup"} (superiority, default) or \code{"ni"} (non-inferiority).
#' @author Maia Muresan (DCR Bern)
#' @importFrom stats approx pnorm qnorm spline uniroot
#' @return An object of class "dx_bin" containing:
#' \itemize{
#'   \item design specifications(design specifications (alpha level, testing framework)
#'   \item required total sample size per group,
#'   \item subgroup sample sizes (diseased / non-diseased),
#'   \item target overall power
#'   \item allocated endpoint-specific powers.
#'
#' }
#' @export
#' @examples
#' unpaired_ss_opt_bin(alpha = 0.05,power = 0.80,se.c = 0.75,se.e = 0.85,sp.c = 0.65,sp.e = 0.75,prev = 0.60)
#' unpaired_ss_opt_bin(alpha = 0.05,power = 0.80,se.c = 0.75,se.e = 0.85,sp.c = 0.65,sp.e = 0.75,prev = 0.60, delta.se=0.05,type.se="ni")
#' unpaired_ss_opt_bin(alpha = 0.05,power = 0.80,se.c = 0.75,se.e = 0.85,sp.c = 0.65,sp.e = 0.75,prev = 0.60, delta.sp=0.05,type.sp="ni")
#'
#'
unpaired_ss_opt_bin <- function (power,
                                  prev,
                                  se.c, se.e,
                                  sp.c,sp.e,
                                  alpha,
                                  delta.se=0,
                                  delta.sp=0,
                                  type.se=c("sup","ni"),
                                  type.sp=c("sup","ni")) {
  # parameter description
  # alpha: desired type I error rate per endpoint
  # power: desired overall power
  # se.c: sensitivity of the comparator test
  # se.e: sensitivity of the experimental test
  # sp.c: specificity of the comparator test
  # sp.e: specificity of the experimental test
  # prev: prevalence



  type.se <- match.arg(type.se)
  type.sp <- match.arg(type.sp)
  # the variance function
  # variance function
  v <- function(c,e){
    c*(1-c) + e*(1-e)
  }


  # ---------- endpoint-level sample size ----------
  n_endpoint <- function(alpha, beta, c, e, delta, type){

    if(type == "sup"){
      z.alpha <- qnorm(alpha/2)
      effect  <- c - e
      V0 <- v(c,c)
      VA <- v(c,e)
      n <- ((z.alpha*sqrt(V0) + qnorm(beta)*sqrt(VA))^2) / (effect^2)
    }

    if(type == "ni"){
      z.alpha <- qnorm(alpha)
      effect  <- c - e - delta
      VA <- v(c,e)
      n <- ((z.alpha + qnorm(beta))^2 * VA) / (effect^2)
    }

    return(ceiling(n))
  }
  f <- function(beta.se){

    beta.sp <- (power + beta.se - 1)/(beta.se - 1)

    if(beta.sp <= 0 || beta.sp >= 1) return(Inf)

    n.se <- n_endpoint(alpha, beta.se,
                       se.c, se.e,
                       delta.se, type.se)

    n.sp <- n_endpoint(alpha, beta.sp,
                       sp.c, sp.e,
                       delta.sp, type.sp)

    N.se <- n.se / prev
    N.sp <- n.sp / (1-prev)

    return(N.se - N.sp)
  }
  beta.se <- uniroot(f,
                     lower = 1e-8,
                     upper = 1-power-1e-8)$root

  beta.sp <- (power + beta.se - 1)/(beta.se - 1)

  # ---------- final sample sizes ----------
  n.se <- n_endpoint(alpha, beta.se,
                     se.c, se.e,
                     delta.se, type.se)

  n.sp <- n_endpoint(alpha, beta.sp,
                     sp.c, sp.e,
                     delta.sp, type.sp)

  N.se <- n.se / prev
  N.sp <- n.sp / (1-prev)

  N.per.group <- ceiling(max(N.se, N.sp))

  # achieved power check
  power.se <- 1 - pnorm(
    (sqrt(N.per.group*prev)*(se.c-se.e-delta.se*(type.se=="ni")) -
       qnorm(ifelse(type.se=="sup",alpha/2,alpha))*sqrt(v(se.c,se.c))) /
      sqrt(v(se.c,se.e))
  )

  power.sp <- 1 - pnorm(
    (sqrt(N.per.group*(1-prev))*(sp.c-sp.e-delta.sp*(type.sp=="ni")) -
       qnorm(ifelse(type.sp=="sup",alpha/2,alpha))*sqrt(v(sp.c,sp.c))) /
      sqrt(v(sp.c,sp.e))
  )

  power.total <- power.se * power.sp

  result <- list(
    alternative="two-sided",
    alpha = alpha,
    hypothesis=list(
      sensitivity = type.se,
      specificity = type.sp
    ),
    design  = "unpaired comparative design",
    analysis = "samplesize",
    target_power = power.total,
    power_se = power.se,
    power_sp = power.sp,
    N_se = n.se,
    N_sp = n.sp,
    N_tot = N.per.group,
    prev = prev
  )

  class(result) <- "dx_bin"
  return(result)
}

