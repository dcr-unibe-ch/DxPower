#' Optimal Sample Size Calculation for an Paired Comparative
#' Diagnostic Accuracy Study with Co-Primary Endpoints
#'
#' @description
#' The function determines the minimal total sample size per group that guarantees the desired overall power for a paired
#' comparative diagnostic accuracy study with co-primary endpoints
#' (sensitivity and specificity) using the intersection–union testing
#' framework (IUT).
#'
#' In a paired design, each subject receives the experimental test, the comparator test and the reference standard.
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
#' \strong{Required assumptions (nuisance parameters)}
#'
#' In addition to sensitivities and specificities, the paired design
#' requires assumptions on:
#'
#' \itemize{
#'   \item Disease prevalence: \eqn{\pi}
#'   \item Discordant proportion among diseased subjects: \eqn{\psi_D}
#'   \item Discordant proportion among non-diseased subjects: \eqn{\psi_{ND}}
#' }
#'
#' These are nuisance parameters and must be specified based on prior
#' studies or clinical knowledge.
#'
#' For diseased subjects:
#'
#' \deqn{
#' |Se_C - Se_E| \le \psi_D \le Se_C + Se_E - 2 Se_C Se_E
#' }
#'
#' For non-diseased subjects:
#'
#' \deqn{
#' |Sp_C - Sp_E| \le \psi_{ND} \le Sp_C + Sp_E - 2 Sp_C Sp_E
#' }
#'
#' Minimal discordance values are therefore:
#'
#' \deqn{\psi_{D,min} = |Se_C - Se_E|}
#' \deqn{\psi_{ND,min} = |Sp_C - Sp_E|}
#'
#' In practice, the sample size may be calculated using these
#' minimal values together with an assumed prevalence \eqn{\pi_{ass}}.
#'
#' \strong{Integration of prevalence}
#'
#' Let \eqn{\theta_{Se,E}} and \eqn{\theta_{Se,C}} denote the sensitivities
#' of the experimental and comparator tests, and
#' \eqn{\theta_{Sp,E}} and \eqn{\theta_{Sp,C}} the corresponding specificities.
#'
#' With total sample size \eqn{N} and disease prevalence \eqn{\pi},
#' subgroup sizes are:
#'
#' \deqn{N_{Se} = N\pi, \qquad N_{Sp} = N(1-\pi).}
#'
#' Sensitivity is evaluated in the diseased subgroup and specificity
#' in the non-diseased subgroup. Assuming independence between these
#' subgroups, overall power under the IUT framework equals:
#'
#' \deqn{Power_{overall} = Power_{Se} \times Power_{Sp}.}
#'
#' The global null hypothesis is the union of the endpoint-specific
#' null hypotheses and is rejected only if both endpoints are significant.
#'
#' To achieve the desired overall power, the function solves numerically:
#'
#' \deqn{(1-\beta_{Se})(1-\beta_{Sp}) = power}
#'
#' under the constraint:
#'
#' \deqn{N_{Se} = N_{Sp},}
#'
#' ensuring optimal allocation of type II error between endpoints.
#'
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
#' @seealso \code{\link{ss_opt_bin}}
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
#' @param psi.d discordant test results for diseased patients
#' @param psi.nd discordant test results for non-diseased patients
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
#'
paired_ss_opt_bin <- function(power,
                              prev,
                              se.c, se.e,
                              sp.c, sp.e,
                              psi.d,
                              psi.nd,
                              alpha,
                              delta.se=0,
                              delta.sp=0,
                              type.se=c("sup","ni"),
                              type.sp=c("sup","ni")) {


  type.se <- match.arg(type.se)
  type.sp <- match.arg(type.sp)

  # -----------------------------
  # Input validation
  # -----------------------------
  if (power <= 0 || power >= 1)stop("power must be between 0 and 1.")
  if (prev <= 0 || prev >= 1)stop("Prevalence must be between 0 and 1.")
  if (se.e <= se.c)stop("se.e must be greater than se.c for superiority testing.")
  if (sp.e <= sp.c)stop("sp.e must be greater than sp.min for superiority testing.")
  if (alpha <= 0 || alpha >= 1)stop("alpha must be between 0 and 1.")

  # sample size calculation for one endpoint following Miettinen (1968)
  n_endpoint <- function(alpha, beta, c, e, psi,  delta=0, type){
    z.beta <- qnorm(1 - beta)
    z.alpha <- qnorm(1-alpha/2)

    if(type == "sup"){
      diff   <- abs(e - c)
      n <- (z.alpha*psi+z.beta*sqrt( (psi^2)-0.25*((c-e)^2)*(3+psi)) )^2 /(psi*(c-e)^2)
    }

    if(type == "ni"){

      p01<-(psi-e+c)/2
      b_u<-delta*(delta+1)*p01
      a_u<-(-(e-c))*(1-delta)-2*(p01+delta)
      p_bar_01<-(-a_u+sqrt(a_u^2-8*b_u))/4
      w_u<-sqrt(2*p01+(e-c)-(e-c)^2)/sqrt(2*p_bar_01-delta-delta^2)

      n <- ((psi -(e-c)^2)*((z.alpha/w_u)+z.beta)^2)/(-delta-(e-c))^2

    }
    return(ceiling(n))
  }

  # calculate power for one endpoint


  # -----------------------------
  # Root function for IUT
  # -----------------------------
  f_root <- function(beta.se){

    beta.sp <- (power + beta.se - 1)/(beta.se - 1)

    if(beta.sp <= 0 || beta.sp >= 1) return(Inf)

    n.se <- n_endpoint(alpha=alpha, beta=beta.se,
                       c=se.c, e=se.e,
                       psi=psi.d,
                       delta=delta.se,
                       type=type.se)

    n.sp <- n_endpoint(alpha=alpha, beta=beta.sp,
                       c=sp.c, e=sp.e,
                       psi=psi.nd,
                       delta=delta.sp,
                       type=type.sp)

    N.se <- n.se / prev
    N.sp <- n.sp / (1 - prev)

    return(N.se - N.sp)
  }

  beta.se <- uniroot(f_root,
                     lower = 1e-8,
                     upper = 1-power-1e-8)$root

  beta.sp <- (power + beta.se - 1)/(beta.se - 1)

  # -----------------------------
  # Final sample sizes
  # -----------------------------
  n.se <- n_endpoint(alpha=alpha, beta=beta.se,
                     c=se.c, e=se.e,
                     psi=psi.d,
                     delta=delta.se,
                     type=type.se)

  n.sp <- n_endpoint(alpha=alpha, beta=beta.sp,
                     c=sp.c, e=sp.e,
                     psi=psi.nd,
                     delta=delta.sp,
                     type=type.sp)

  N.se <- n.se / prev
  N.sp <- n.sp / (1-prev)

  N<- ceiling(max(N.se, N.sp))

  # -----------------------------
  # Power check
  # -----------------------------
  power_endpoint <- function(n, alpha, c, e, psi, delta=0, type) {

    if(type == "sup"){
      diff <- abs(c-e)
      z <- (sqrt(n*psi)*(diff) - qnorm(1-alpha/2)*psi)/sqrt((psi^2)-0.25*(diff^2)*(3+psi))
      power <- pnorm(z)
    }
    if(type == "ni"){
      z.alpha <- qnorm(1-alpha/2)
      p01<-(psi-e+c)/2
      bu<-delta*(delta+1)*p01
      au<-(e-c)*(1-delta)-2*(p01+delta)
      p_bar_01<-(-au+sqrt(au^2-8*bu))/4
      w_u<-sqrt(2*p01+(e-c)-(e-c)^2)/sqrt(2*p_bar_01-delta-delta^2)
      z <- (-z.alpha/w_u)-(delta+(e-c))/sqrt((psi-(e-c)^2)/n)
      power <- pnorm(z)
    }
    return(power)
  }


  power.se <- power_endpoint(n=N*prev,
                             alpha=alpha,
                             c=se.c, e=se.e,
                             psi=psi.d,
                             delta=delta.se,
                             type=type.se)

  power.sp <- power_endpoint(n=N*(1-prev),
                             alpha=alpha,
                             c=sp.c, e=sp.e,
                             psi=psi.nd,
                             delta=delta.sp,
                             type=type.sp)

  power.total <- power.se * power.sp


  result <- list(
    alternative="two-sided",
    alpha = alpha,
    hypothesis=list(
      sensitivity = type.se,
      specificity = type.sp
    ),
    design  = "paired comparative design",
    analysis = "samplesize",
    target_power = power.total,
    power_se = power.se,
    power_sp = power.sp,
    N_se = n.se,
    N_sp = n.sp,
    N_tot = N,
    prev = prev
  )

  class(result) <- "dx_bin"
  return(result)
}


