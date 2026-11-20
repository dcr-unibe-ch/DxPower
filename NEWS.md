# DxPower 0.1.0

* Initial dev version.

## Functions

* `power.copr.test()`: power or sample size of a single-test design (superiority
  against fixed minimum sensitivity and specificity), with Wald or logit
  approximation. Give exactly one of `N` or `power`.
* `power.copr.comp.test()`: comparative designs, unpaired and paired, with
  superiority or non-inferiority per endpoint. The paired design uses the
  discordance (`psi.d`, `psi.nd`) or the phi correlation (`rho.d`, `rho.nd`) of
  the two tests.
* `paired_bounds()`, `rho_to_psi()` and `psi_to_rho()`: feasible range of the
  discordance and correlation and conversion between them.
* `ss.copr.conv()`: conventional single-test sample size (per-endpoint sample
  size divided by prevalence, maximum over the endpoints).
* `power.copr.emp()` and `power.copr.comp.emp()`: Monte Carlo power for the
  single-test and the unpaired design.

## Exact power

* `exact = TRUE` computes exact binomial power (single-test and unpaired
  designs). This is an **experimental** feature and it
  needs further development and testing.

## Behaviour

* Overall power is the product of the endpoint powers (the diseased and
  non-diseased samples are independent).
* Inputs are validated (probabilities, sample sizes, `n_sim`). A warning is
  given if the effective one-sided level is above 0.025.
* Results are of class `power.htest`.
