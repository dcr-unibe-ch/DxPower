#' Power Calculation for an Unpaired Comparative
#' Diagnostic Accuracy Study with Co-Primary Endpoints

#' @description
#' The function determines the overall power for an unpaired comparative
#' diagnostic accuracy study with co-primary endpoints (sensitivity
#' and specificity) under the intersection–union testing (IUT) framework.
#'
#' Superiority and non-inferiority hypotheses may be specified
#' independently for sensitivity and specificity.
#'
#' Sensitivity is evaluated in the diseased subgroup and specificity in
#' the non-diseased subgroup. Power calculations are based on large-sample
#' normal approximations for the difference between two independent
#' binomial proportions.
#'
#' @seealso \code{\link{unpaired_ss_opt_bin}}
#'
#'
#' @references
#' Stark, M., Hesse, M., Brannath, W., & Zapf, A. (2022).
#' Blinded sample size re-estimation in a comparative diagnostic
#' accuracy study.
#' \emph{BMC Medical Research Methodology}, 22, 115.
#' \doi{10.1186/s12874-022-01564-2}
#'
#'
#' @param N Total sample size.
#' @param prev Disease prevalence (0 < prev < 1).
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
#'
#' @importFrom stats approx pnorm qnorm spline uniroot
#' @return An object of class `"dx_bin"` containing:
#' \itemize{
#'   \item design specifications (design specifications (alpha level, testing framework, CI method)
#'   \item achieved overall power,
#'   \item endpoint-specific powers,
#'   \item subgroup sample sizes.

#' }
#'
#' @export
#' @examples
#' unpaired_power_opt_bin(
#' N=910,
#' alpha=0.05,
#' se.c=0.75,
#' se.e=0.85,
#' sp.c=0.65,
#' sp.e=0.75,
#' prev=0.6
#' )
#' unpaired_power_opt_bin(
#' N=910,
#' alpha=0.05,
#' se.c=0.75,
#' se.e=0.85,
#' sp.c=0.65,
#' sp.e=0.75,
#' prev=0.6,
#' delta.sp=0.05,
#' type.sp="ni"
#' )
#'




unpaired_power_opt_bin <- function(
    N,                 # total sample size per group
    prev,
    alpha,
    se.c, se.e,
    sp.c,sp.e,
    delta.se=0,
    delta.sp=0,
    type.se=c("sup","ni"),
    type.sp=c("sup","ni")) {

  type.se <- match.arg(type.se)
  type.sp <- match.arg(type.sp)

  # the variance function
  v <- function(c,e){
    c*(1-c) + e*(1-e)
  }

  # power for one endpoint
  # ---------- endpoint power ----------
  power_endpoint <- function(n, alpha, c, e, delta, type){
    z.alpha <- qnorm(1-alpha/2)

    if(type == "sup"){
      effect  <- c - e
      V0 <- v(c,c)
      VA <- v(c,e)
      z <- (sqrt(n)*effect - z.alpha*sqrt(V0)) / sqrt(VA)
    }

    if(type == "ni"){
      effect  <- c - e - delta
      VA <- v(c,e)
      z <- ((sqrt(n)*(effect)/sqrt(VA))-z.alpha)
    }

    power <- 1 - pnorm(z)

    return(power)
  }

  # sample sizes per endpoint per arm
  Nd <- N * prev #diseased
  Nh <- N * (1 - prev) #non-diseased(healty)

  # compute powers
  power.se <- power_endpoint(n=Nd, alpha=alpha,
                             c=se.c, e=se.e,
                             delta=delta.se, type=type.se)

  power.sp <- power_endpoint(n=Nh, alpha=alpha,
                             c=sp.c, e=sp.e,
                             delta=delta.sp, type=type.sp)

  power.total <- power.se * power.sp

  result <- list(
    alpha = alpha,
    hypothesis=list(
      sensitivity = type.se,
      specificity = type.sp
    ),
    design  = "unpaired comparative design",
    analysis = "power",
    overall_power = power.total,
    power_se = power.se,
    power_sp = power.sp,
    N_se = Nd,
    N_sp = Nh,
    N_tot = N,
    prev = prev
  )

  class(result) <- "dx_bin"
  return(result)
}


