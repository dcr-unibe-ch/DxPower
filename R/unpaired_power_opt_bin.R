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
#' @param rho Correlation between the sensitivity and specificity estimators
#'   (default 0, i.e. independent endpoints). Overall power is the
#'   probability that both test statistics exceed their critical values under
#'   a bivariate normal distribution with this correlation.
#' @param variance Variance of the difference estimator under the null used
#'   for the critical value: \code{"null"} (default; variance at the comparator
#'   proportion for superiority, alternative variance for non-inferiority),
#'   \code{"pooled"} (pooled proportion \eqn{(c+e)/2}) or \code{"unpooled"}
#'   (alternative variance). The alternative variance is always unpooled.
#' @author Maia Muresan (DCR Bern)
#'
#' @importFrom stats approx pnorm qnorm spline uniroot
#' @importFrom mvtnorm pmvnorm TVPACK
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




unpaired_power_opt_bin <- function(
    N,                 # total sample size per group
    prev,
    alpha,
    se.c, se.e,
    sp.c,sp.e,
    delta.se=0,
    delta.sp=0,
    type.se=c("sup","ni"),
    type.sp=c("sup","ni"),
    rho = 0,
    variance = c("null", "pooled", "unpooled")) {

  type.se <- match.arg(type.se)
  type.sp <- match.arg(type.sp)
  variance <- match.arg(variance)

  if (rho <= -1 || rho >= 1) stop("rho must be strictly between -1 and 1.")

  # the variance function
  v <- function(c,e){
    c*(1-c) + e*(1-e)
  }

  # per-endpoint ingredients of the test "estimated difference > threshold":
  #   mean      = true difference + NI margin (0 for superiority)
  #   var       = variance of the difference estimator under the alternative / n
  #   threshold = z_alpha * sd of the estimator under the null / sqrt(n)
  endpoint_params <- function(n, alpha, c, e, delta, type){
    z.alpha <- qnorm(1-alpha/2)
    margin  <- if(type == "ni") delta else 0

    VA <- v(c,e)
    V0 <- switch(variance,
      pooled   = {p.bar <- (c + e)/2; 2*p.bar*(1-p.bar)},
      unpooled = VA,
      null     = if(type == "sup") v(c,c) else VA)

    list(mean = (e - c) + margin,
         var  = VA / n,
         thr  = z.alpha * sqrt(V0 / n))
  }

  # sample sizes per endpoint per arm
  Nd <- N * prev #diseased
  Nh <- N * (1 - prev) #non-diseased(healty)

  par.se <- endpoint_params(n=Nd, alpha=alpha,
                            c=se.c, e=se.e,
                            delta=delta.se, type=type.se)
  par.sp <- endpoint_params(n=Nh, alpha=alpha,
                            c=sp.c, e=sp.e,
                            delta=delta.sp, type=type.sp)

  # marginal powers
  power.se <- 1 - pnorm(par.se$thr, par.se$mean, sqrt(par.se$var))
  power.sp <- 1 - pnorm(par.sp$thr, par.sp$mean, sqrt(par.sp$var))

  # joint power: P(both estimated differences exceed their thresholds)
  # under a bivariate normal with diagonal binomial variances and
  # correlation rho between the sensitivity and specificity estimators
  K <- 2
  R <- diag(K)
  R[1, 2] <- R[2, 1] <- rho

  if (min(eigen(R, symmetric = TRUE)$values) <= 0)
    stop("R is not positive definite; choose different rho's.")

  V <- diag(c(par.se$var, par.sp$var))
  D <- diag(sqrt(diag(V)))
  Sigma <- D %*% R %*% D

  power.total <- as.numeric(
    mvtnorm::pmvnorm(lower = c(par.se$thr, par.sp$thr),
                     upper = c(Inf, Inf),
                     mean  = c(par.se$mean, par.sp$mean),
                     sigma = Sigma,
                     algorithm = mvtnorm::TVPACK())
  )

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
    rho = rho,
    N_se = Nd,
    N_sp = Nh,
    N_tot = N,
    prev = prev
  )

  class(result) <- "dx_bin"
  return(result)
}
