# Smallest integer N with pwr_fun(N) >= power, for a power function that is
# non-decreasing in N: bracket by doubling, then bisect.
find_n <- function(pwr_fun, power, n_low = 2, n_high = 1000, n_max = 1e7) {
  if (pwr_fun(n_low) >= power) {
    warning("Target power is already reached at n_low = ", n_low, ".")
    return(n_low)
  }

  while (pwr_fun(n_high) < power) {
    n_low <- n_high
    n_high <- n_high * 2
    if (n_high > n_max) stop("Upper bound exploded; check inputs.")
  }

  lo <- n_low
  hi <- n_high
  while (lo + 1 < hi) {
    mid <- floor((lo + hi) / 2)
    if (pwr_fun(mid) >= power) hi <- mid else lo <- mid
  }
  hi
}
