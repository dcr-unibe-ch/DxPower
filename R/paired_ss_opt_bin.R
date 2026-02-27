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
#' @examples
#' paired_ss_opt_bin(alpha=0.05, power=0.80, se.c=0.81, se.e=0.90,sp.c=0.66, sp.e=0.80,psi.d=0.09, psi.nd=0.14,prev=0.05)
#' paired_ss_opt_bin(alpha=0.05, power=0.80, se.c=0.81, se.e=0.90,sp.c=0.66, sp.e=0.80,psi.d=0.09, psi.nd=0.14,prev=0.05, delta.se=0,delta.sp=0, type.se="sup",type.sp="sup")
#' paired_ss_opt_bin(alpha=0.05, power=0.80, se.c=0.81, se.e=0.90,sp.c=0.66, sp.e=0.80,psi.d=0.09, psi.nd=0.20, prev=0.5, delta.sp = 0.10,type.sp = "ni")
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
  # ==========================================================
  # Co-primary sample size for paired diagnostic study
  # Framework:
  # -paired comparison (same subjects receive boths test in adidtion to the ref.standard)
  # -co-primary endpoints: Sens and Spec
  # -Intersection-Union Test(IUT) => Overall power=Power_Se x Power_Sp
  # IUT hypothesis:
  # - Superioroty: H0_global: d_se<=0 U d_sp<=0 (vs H1_global: d_se>0 AND d_sp>0)
  #   *sample size via Maria stark (2022) asymptotic formula
  # - Non-inferiority: H0_global: d_sp<=0 U d_se<=-Δ (vs H1_gloval: d_se>0 AND d_se > -Δ)
  #   *Sample size via Liu et al. score-based derivation
  # where d=e-c (experimental - comparator test)
  #       ψ = discordance probability = p01 + p10
  # ==========================================================

  type.se <- match.arg(type.se)
  type.sp <- match.arg(type.sp)

  # -----------------------------
  # Input validation
  # -----------------------------
  if (power <= 0 || power >= 1)stop("power must be between 0 and 1.")
  if (prev <= 0 || prev >= 1)stop("Prevalence must be between 0 and 1.")
  if (alpha <= 0 || alpha >= 1)stop("alpha must be between 0 and 1.")
  # -----------------------------
  # Hypothesis consistency checks
  # -----------------------------
  # Superiority requires expected improvement
  if(type.se=="sup" && se.e <= se.c)stop("Sensitivity superiority requires se.e > se.c")
  if(type.sp=="sup" && sp.e <= sp.c)stop("Specificity superiority requires sp.e > sp.c.")
  # Non-inferiority requires true effect above NI boundary
  if(type.se=="ni" && (se.e - se.c) <= -delta.se)
    stop("Sensitivity NI invalid: true difference (se.e - se.c) must be > -delta.se")
  if(type.sp=="ni" && (sp.e - sp.c) <= -delta.sp)
    stop("Specificity NI invalid: true difference (sp.e - sp.c) must be > -delta.sp")



  # Informative warning: large NI margin makes test trivial

  if(type.se=="ni" && delta.se > abs(se.e - se.c)){
    warning(
      paste0(
        "Sensitivity NI margin (delta.sp = ", delta.se,
        ") is larger than the expected difference (sp.e - sp.c = ", round(se.e - se.c,3), "). ",
        "NI objective is far from boundary and unlikely to drive sample size."
      )
    )
  }
  if(type.sp=="ni" && delta.sp > abs(sp.e - sp.c)){
    warning(
      paste0(
        "Specificity NI margin (delta.sp = ", delta.sp,
        ") is larger than the expected difference (sp.e - sp.c = ",
        round(sp.e - sp.c, 3), "). ",
        "NI objective is far from boundary and unlikely to drive sample size."
      )
    )
  }

  # -----------------------------
  # Discordance bounds check
  # ψ must satisfy theoretical bounds:
  # |d| ≤ ψ ≤ c + e − 2ce
  # -----------------------------
  # Sensitivity bounds
  psi.d.lower <- abs(se.c - se.e)
  psi.d.upper <- se.c + se.e - 2*se.c*se.e

  if(psi.d < psi.d.lower || psi.d > psi.d.upper){
    stop(
      paste0("psi.d outside theoretical bounds for paired binary data.",
        " psi.d must satisfy: ",
        round(psi.d.lower, 4), " <= psi.d <= ",
        round(psi.d.upper, 4),
        ". You provided psi.d = ", psi.d
      )
    )
  }

  # Specificity bounds
  psi.nd.lower <- abs(sp.c - sp.e)
  psi.nd.upper <- sp.c + sp.e - 2*sp.c*sp.e

  if(psi.nd < psi.nd.lower || psi.nd > psi.nd.upper){
    stop(
      paste0("psi.nd outside theoretical bounds for paired binary data.",
        "psi.nd must satisfy: ",
        round(psi.nd.lower, 4), " <= psi.nd <= ",
        round(psi.nd.upper, 4),
        ". You provided psi.nd = ", psi.nd
      )
    )
  }

  # ==========================================================
  # Sample size for a single endpoint
  # ==========================================================
  n_endpoint <- function(alpha, beta, c, e, psi,  delta=0, type){

    z.beta <- qnorm(1 - beta)
    z.alpha <- qnorm(1-alpha/2)
    d   <- e - c # difference in paired proportions

    # Superiority: Maria stark (2022)
    # -----------------------------
    if(type == "sup"){
      n <- (z.alpha*psi+z.beta*sqrt((psi^2)-0.25*((d^2)*(3+psi))))^2 /(psi*(d^2))
    }

    # Non-inferiority: Liu score method
    # Requires quadratic solution for p_bar_01 (p̄01)
    # -----------------------------
    if(type == "ni"){
      # ---- NEW CHECK 1: Alternative must hold ----
      if(delta + d <= 0)
        stop("NI alternative not satisfied: delta + (e - c) must be > 0.")

      p01<-(psi-d)/2
      b_u<-delta*(delta+1)*p01
      a_u<-(-d)*(1-delta)-2*(p01+delta)
      # ---- NEW CHECK 2: Discriminant positive ----
      disc <- a_u^2 - 8*b_u
      if(disc <= 0)
        stop(paste0(
          "Invalid NI configuration: no real solution exists for the assumed ",
          "margin (delta = ", delta,
          ") and discordance (psi = ", psi, "). ",
          "The margin may be too small or psi too close to its boundary."
        ))

      p_bar_01<-(-a_u+sqrt(disc))/4
      w_u<-sqrt(2*p01+d-d^2)/sqrt(2*p_bar_01-delta-delta^2)


      n <- ((psi -d^2)*(z.beta+z.alpha/w_u)^2)/(delta+d)^2

    }
    return(n)
  }

  # ==========================================================
  # Beta allocation under IUT:
  # (1 - beta.se)(1 - beta.sp) = overall.power
  # ==========================================================
  beta.grid <- seq(1e-6, 1 - power - 1e-6, length.out = 10000)

  diff.values <- vapply(beta.grid, function(beta.se){

    beta.sp <- 1 - power/(1 - beta.se)

    if(beta.sp <= 0 || beta.sp >= 1) return(Inf)

    n.se <- n_endpoint(alpha, beta.se,
                       se.c, se.e,
                       psi.d,
                       delta.se,
                       type.se)

    n.sp <- n_endpoint(alpha, beta.sp,
                       sp.c, sp.e,
                       psi.nd,
                       delta.sp,
                       type.sp)

    N.se <- n.se / prev
    N.sp <- n.sp / (1 - prev)

    abs(N.se - N.sp)

  }, numeric(1))

  beta.se <- beta.grid[which.min(diff.values)]
  beta.sp <- 1-power/(1-beta.se)

  # -----------------------------
  # Final sample sizes
  # -----------------------------
  #sample size of disease individuals
  n.se <- n_endpoint(alpha=alpha, beta=beta.se,
                     c=se.c, e=se.e,
                     psi=psi.d,
                     delta=delta.se,
                     type=type.se)
  #sample size of non-disease individuals
  n.sp <- n_endpoint(alpha=alpha, beta=beta.sp,
                     c=sp.c, e=sp.e,
                     psi=psi.nd,
                     delta=delta.sp,
                     type=type.sp)

  N.se <- n.se / prev #total sample size including at least n_se individuals
  N.sp <- n.sp / (1-prev)#total sample size including at least n_sp individuals

  N<- max(N.se, N.sp)

  # ==========================================================
  # Power calculation (consistency check)
  # ==========================================================
  power_endpoint <- function(n, alpha, c, e, psi, delta=0, type) {

    z.alpha<-qnorm(1-alpha/2)
    d<-e-c

    if(type == "sup"){
      z <- (sqrt(n*psi)*(d) - z.alpha*psi)/sqrt((psi^2)-0.25*(d^2)*(3+psi))
      power <- pnorm(z)
    }
    if(type == "ni"){
      p01<-(psi-d)/2
      b_u<-delta*(delta+1)*p01
      a_u<-(-d)*(1-delta)-2*(p01+delta)
      p_bar_01<-(-a_u+sqrt(a_u^2-8*b_u))/4
      w_u<-sqrt(2*p01+d-d^2)/
        sqrt(2*p_bar_01-delta-delta^2)
      sigma<-sqrt((psi-d^2)/n)
      z <- (z.alpha/w_u)-(delta+d)/sigma
      power <- 1-pnorm(z)
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






