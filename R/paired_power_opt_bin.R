#' Power Calculation for a Paired Comparative
#' Diagnostic Accuracy Study with Co-Primary Endpoints
#'
#' @description
#' The function determines the achieved overall power for a fixed total sample size \eqn{N} in a paired
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
#' Sensitivity is evaluated in the diseased subgroup and specificity
#' in the non-diseased subgroup.
#'
#' Let \eqn{\theta_{Se,E}} and \eqn{\theta_{Se,C}} denote the sensitivities
#' of the experimental and comparator tests, and
#' \eqn{\theta_{Sp,E}} and \eqn{\theta_{Sp,C}} the corresponding specificities.
#' With disease prevalence \eqn{\pi}, subgroup sizes are:
#'
#' \deqn{N_{Se} = N \pi}
#' \deqn{N_{Sp} = N (1-\pi)}
#'
#' \strong{Intersection–Union Testing (IUT)}
#'
#' The global null hypothesis is the union of the endpoint-specific
#' null hypotheses. Study success requires rejection of both
#' endpoint-specific null hypotheses.
#'
#' Under independence of diseased and non-diseased subgroups,
#' the overall power is:
#'
#' \deqn{Power_{overall} = Power_{Se} \times Power_{Sp}.}
#'
#' No multiplicity adjustment of \eqn{\alpha} is required under the
#' IUT framework.
#'
#'
#' \strong{Endpoint-specific sample size formulas}
#'
#'  Let \eqn{d = \theta_E - \theta_C} denote the difference between the
#' experimental and comparator test (sensitivity or specificity,
#' depending on the endpoint).
#'
#' For each endpoint, power is calculated from the
#' large-sample normal distribution of the corresponding
#' paired test statistic.
#'
#' For superiority hypotheses (\eqn{H_0: d \le 0}),
#' power is calculated using the large-sample normal
#' approximation for paired binary data based on the
#' asymptotic variance expression described by
#' Miettinen (1968) and applied in Stark et al. (2022).
#'
#' For non-inferiority hypotheses, power is calculated using the
#' score test for paired binary data described by Liu et al. (2002).
#' The test is based on a restricted maximum likelihood (RMLE)
#' estimator under the non-inferiority null hypothesis
#' \eqn{H_0: d = -\Delta}, so that the variance of the test statistic
#' is evaluated under this null constraint. Power is then obtained using the large-sample normal
#' approximation under the assumed true effect (alternative)
#' \eqn{H_1: d > -\Delta}.
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
#' Discordant parameters must satisfy the theoretical bounds.
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
#'
#' \strong{Important assumptions}
#'
#' \itemize{
#'   \item Asymptotic normal approximation is used.
#'   \item Diseased and non-diseased subgroups are treated as independent.
#'   \item The true effect must satisfy the alternative hypothesis
#'         (e.g., \eqn{d + \Delta > 0} for non-inferiority).
#' }
#'
#'
#' The function evaluates the achieved endpoint-specific powers for a fixed \eqn{N},
#' and multiplies them according to the IUT framework.
#' No numerical optimisation or beta allocation is performed.
#'
#'
#' @seealso
#' \code{\link{paired_ss_opt_bin}} for sample size determination,
#' \code{\link{ss_opt_bin}} and
#' \code{\link{power_opt_bin}} for single-test design
#' \code{\link{unpaired_ss_opt_bin}} and
#' \code{\link{unpaired_power_opt_bin}} for unpaired designs.
#' @references
#' Stark, M., Hesse, M., Brannath, W., & Zapf, A. (2022).
#' Blinded sample size re-estimation in a comparative diagnostic
#' accuracy study.
#' \emph{BMC Medical Research Methodology}, 22, 115.
#' \doi{10.1186/s12874-022-01564-2}
#'
#' Liu, J. P., Hsueh, H. M., Hsieh, E., & Chen, J. J. (2002).
#' Tests for equivalence or non-inferiority for paired binary data.
#' \emph{Statistics in Medicine}, 21(2), 231–245.
#' \doi{10.1002/sim.1012}
#'
#' Miettinen, O.S. (1968).
#' The matched pairs design in the case of all-or-none responses.
#' \emph{Biometrics}, 24(2), 339–352.
#' \doi{10.2307/2528039}
#'
#' @param N Total sample size.
#' @param prev Disease prevalence (0 < prev < 1).
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
#' @importFrom stats  pnorm qnorm
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
#'
#' # Example 1: Superiority for both sensitivity and specificity
#' # Compute achieved power for a fixed total sample size
#' paired_power_opt_bin(
#'   N = 134,            # total sample size
#'   prev = 0.47,        # disease prevalence
#'   se.c = 0.81,        # comparator sensitivity
#'   se.e = 0.90,        # experimental sensitivity
#'   sp.c = 0.66,        # comparator specificity
#'   sp.e = 0.80,        # experimental specificity
#'   psi.d = 0.09,       # discordance (diseased)
#'   psi.nd = 0.14,      # discordance (non-diseased)
#'   alpha = 0.05
#' )
#'
#' # Explicit superiority specification (same as default)
#' paired_power_opt_bin(
#'   N = 134,
#'   prev = 0.47,
#'   se.c = 0.81, se.e = 0.90,
#'   sp.c = 0.66, sp.e = 0.80,
#'   psi.d = 0.09,
#'   psi.nd = 0.14,
#'   alpha = 0.05,
#'   type.se = "sup",
#'   type.sp = "sup"
#' )
#'
#'
#' # Example 2: Sensitivity superiority and specificity non-inferiority
#' # Non-inferiority margin for specificity = 0.10
#' paired_power_opt_bin(
#'   N = 134,
#'   prev = 0.47,
#'   se.c = 0.81, se.e = 0.90,
#'   sp.c = 0.66, sp.e = 0.80,
#'   psi.d = 0.09,
#'   psi.nd = 0.14,
#'   alpha = 0.05,
#'   delta.sp = 0.10, # non-inferiority margin for specificity
#'   type.sp = "ni"
#' )

paired_power_opt_bin <- function(N,
                                prev,
                                se.c, se.e,
                                sp.c, sp.e,
                                psi.d,
                                psi.nd,
                                alpha,
                                delta.se = 0,
                                delta.sp = 0,
                                type.se = c("sup","ni"),
                                type.sp = c("sup","ni")) {

  # ==========================================================
  # Power calculation for paired co-primary diagnostic study
  # Fixed total sample size N
  #
  # Overall power (IUT):
  #   Power_total = Power_Se × Power_Sp
  #
  # n_Se = N × prevalence
  # n_Sp = N × (1 - prevalence)
  # ==========================================================

  type.se <- match.arg(type.se)
  type.sp <- match.arg(type.sp)

  # -----------------------------
  # Input validation
  # -----------------------------
  if (N <= 0) stop("Total sample size N must be positive.")
  if (prev <= 0 || prev >= 1) stop("Prevalence must be between 0 and 1.")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1.")
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
        "Non-inferiority may be easy to demonstrate with this margin."
      )
    )
  }
  if(type.sp=="ni" && delta.sp > abs(sp.e - sp.c)){
    warning(
      paste0(
        "Specificity NI margin (delta.sp = ", delta.sp,
        ") is larger than the expected difference (sp.e - sp.c = ",
        round(sp.e - sp.c, 3), "). ",
        "Non-inferiority may be easy to demonstrate with this margin."
      )
    )
  }

  d.se <- se.e - se.c
  d.sp <- sp.e - sp.c

  # sample sizes in each stratum
  n.se <- N * prev
  n.sp <- N * (1 - prev)

  # ==========================================================
  # Endpoint-specific power
  # ==========================================================
  power_endpoint <- function(n, alpha, c, e, psi, delta=0, type){

    z.alpha <- qnorm(1 - alpha/2)
    d <- e - c

    # -----------------------------
    # Superiority (Maria Strak (2022))
    # -----------------------------
    if(type == "sup"){

      z <- (sqrt(n*psi)*d - z.alpha*psi) /
        sqrt(psi^2 - 0.25*d^2*(3+psi))
      power<-pnorm(z)

    }

    # -----------------------------
    # Non-inferiority (Liu score)
    # -----------------------------
    if(type == "ni"){

      # ---- NEW CHECK 1: Alternative must hold ----
      if(delta + d <= 0)
        stop("NI alternative not satisfied: delta + (e - c) must be > 0.")


      p01 <- (psi - d)/2
      b_u <- delta*(delta+1)*p01
      a_u <- (-d)*(1-delta) - 2*(p01+delta)

      disc <- a_u^2 - 8*b_u
      if(disc <= 0)
        stop(paste0(
          "Invalid NI configuration: no real solution exists for the assumed ",
          "margin (delta = ", delta,
          ") and discordance (psi = ", psi, "). ",
          "The margin may be too small or psi too close to its boundary."
        ))

      p_bar_01 <- (-a_u + sqrt(disc))/4

      w_u <- sqrt(2*p01 + d - d^2) /
        sqrt(2*p_bar_01 - delta - delta^2)

      sigma <- sqrt((psi - d^2)/n)

      z <- (z.alpha/w_u) - (delta + d)/sigma
      power<-1-pnorm(z)

    }
    return(power)
  }

  power.se <- power_endpoint(n=n.se, alpha=alpha,
                             c=se.c, e=se.e,
                             psi=psi.d,
                             delta=delta.se,
                             type=type.se)

  power.sp <- power_endpoint(n=n.sp, alpha=alpha,
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
    analysis = "power",
    overall_power = power.total,
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

