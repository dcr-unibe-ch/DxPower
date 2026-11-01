
<!-- README.md is generated from README.Rmd. Please edit that file -->

``` r
devtools::load_all()
#> ℹ Loading DxPower
```

# DxPower

<!-- badges: start -->
<!-- badges: end -->

The goal of DxPower is to provide functions for power and sample size
calculation in confirmatory diagnostic accuracy trials.

The methodology is based on asymptotic normal approximations and an
intersection–union testing framework, where both co-primary endpoints
(i.e. sensitivity and specificity) must reject their respective null
hypotheses.

## Contents

- [Installation](#installation)
- [Current Implementation Status](#current-implementation-status)
- [Statistical Framework](#statistical-framework)
  - [Study Setting](#study-setting)
  - [Hypotheses](#hypotheses)
  - [Test Statistic and
    Approximation](#test-statistic-and-approximation)
  - [Overall Power Definition](#overall-power-definition)
- [Functions](#functions)
  - [Overview](#overview)
  - [Sample size calculation](#sample-size)
  - [Power calculation](#power-calculation)
  - [Comparative designs](#comparative-designs)
  - [Conventional sample size](#conventional)
  - [Empirical (simulated) power](#empirical)
- [Example](#example)
- [References](#references)

## Installation

You can install the development version of DxPower from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("CTU-Bern/DxPower")
```

## Current Implementation Status

The current version implements sample size determination and power
evaluation for confirmatory diagnostic accuracy studies with:

- Binary endpoints (sensitivity and specificity)

- Single-test designs (superiority against fixed null values) and
  comparative designs (unpaired and paired, superiority or
  non-inferiority)

- Co-primary testing using an intersection–union framework

- Wald-type (or logit) asymptotic normal approximations, and exact
  binomial power for the single-test and the unpaired comparative design

Extensions to continuous endpoints and more general co-primary
objectives are planned.

## Statistical Framework

### Study Setting

We consider a confirmatory diagnostic accuracy study with:

- Sensitivity $Se$
- Specificity $Sp$
- Disease prevalence $\pi$
- Total sample size $N$

The number of diseased (cases) and non-diseased (control) subjects is:

$$N_{Se}=N\pi,  N_{Sp}=N(1-\pi)$$

Diseased and non-diseased subjects form disjoint subgroups. Under
standard binomial sampling assumptions, the corresponding estimators of
sensitivity and specificity are asymptotically independent.

### Hypotheses

We consider one-sided superiority hypotheses:

For sensitivity:

$$
H_{0,Se}: Se \le Se_0
\quad\text{vs}\quad
H_{1,Se}: Se > Se_0
$$

For specificity:

$$
H_{0,Sp}: Sp \le Sp_0
\quad\text{vs}\quad
H_{1,Sp}: Sp > Sp_0
$$

Each endpoint is tested at level $\alpha = 0.025$.

Because both null hypotheses must be rejected to declare overall
success, an intersection–union test (IUT) framework is used. Under the
global null hypothesis the overall Type I error rate is controlled at
level $\alpha$ without multiplicity adjustment.

### Test Statistic and Approximation

Both endpoints are evaluated using a Wald-type normal approximation.
Under large-sample theory:

$$\hat{Se} \sim N\left( Se,\; \frac{Se(1-Se)}{N_{Se}} \right),\hat{Sp} \sim N\left( Sp,\; \frac{Sp(1-Sp)}{N_{Sp}} \right) $$

The sample size formula uses the classical normal approximation with:

- Variance under the null:
  $Var(Se_0)=Se_0(1-Se_0), Var(Sp_0)=Sp_0(1-Sp_0)$
- Variance under the alternative:
  $Var(Se_1)=Se_1(1-Se_1), Var(Sp_1)=Sp_1(1-Sp_1)$

### Overall Power Definition

Because both endpoints must reject the null hypothesis, an
intersection–union test (IUT) is applied.

Assuming asymptotic independence:

$$
Power_{Se} \times Power_{Sp}=Power_{overall} 
$$

$$
(1-\beta_{Se}) \times (1-\beta_{Sp})=(1-\beta_{Overall})
$$

Individual endpoint powers are therefore typically larger than the
target overall power.

## Functions

### Overview

| Function                 | Purpose                                                              |
|--------------------------|----------------------------------------------------------------------|
| `power.copr.test()`      | Sample size (`power` given) or power (`N` given), single-test design |
| `power.copr.comp.test()` | Sample size or power, comparative unpaired or paired design          |
| `ss.copr.conv()`         | Conventional per-endpoint sample size, single-test design            |
| `power.copr.emp()`       | Monte Carlo power, single-test design                                |
| `power.copr.comp.emp()`  | Monte Carlo power, unpaired comparative design                       |

All sample size and power functions return an object of class
`power.htest`.

### Sample size calculation

`power.copr.test(power = ...)` computes the required total sample size
for a desired overall power and true prevalence $\pi$ in a confirmatory
diagnostic accuracy study with co-primary endpoints.

#### Sample Size Formula

For sensitivity:

$$
N_{Se} =
\frac{\left(
z_{1-\alpha} \sqrt{Se_0(1-Se_0)}
+
z_{1-\beta_{Se}} \sqrt{Se_1(1-Se_1)}
\right)^2}
{\pi \times (Se_1 - Se_0)^2}
$$

For specificity:

$$
N_{Sp} =
\frac{\left(
z_{1-\alpha} \sqrt{Sp_0(1-Sp_0)}
+
z_{1-\beta_{Sp}} \sqrt{Sp_1(1-Sp_1)}
\right)^2}
{(1-\pi) \times (Sp_1 - Sp_0)^2}
$$

where:

- $z_{1-\alpha}$ is the one-sided critical value  
- $z_{1-\beta}$ corresponds to the endpoint-specific power  
- $Se_0, Sp_0$ are null values  
- $Se_1, Sp_1$ are alternative values
- $\pi$ is the prevalence; determines the allocation of diseased and
  non-diseased subjects.

#### Root-Finding Condition

The study must satisfy two conditions simultaneously:

1.  The overall power constraint

$$
(1-\beta_{Se}) \times (1-\beta_{Sp})=(1-\beta_{Overall})
$$

2.  Equal total sample size requirement $$
    N_{Se} = N_{Sp}
    $$

Since the two endpoint-specific powers are linked through the overall
power constraint, only one free parameter remains (e.g.,$\beta_{Se}$).

$$
\beta_{Se} = 1- \frac{1-\beta_{Overall}}{1-\beta_{Sp}}
$$

The algorithm therefore:

- Express $\beta_{Sp}$ as a function of $\beta_{Se}$
- Substitutes into the two sample size formulas
- Numerically solves for the value of $\beta_{Se}$ such that

$$
N_{Se}(\beta_{Se}) = N_{Sp}(\beta_{Sp})
$$

This ensures that:

- Both endpoints require the same total sample size
- The overall power equals the desired target
- The study is optimally balanced

The resulting total sample size is $N_{total} = N_{Se} + N_{Sp}$ with
allocation driven by prevalence.

### Power calculation

`power.copr.test(N = ...)` computes the achieved overall power for a
fixed total sample size and true prevalence $\pi$ in a confirmatory
diagnostic accuracy study with co-primary endpoints.

#### Endpoint - Specific Power Formula

Given total $N$ and true prevalence $\pi$:

$$
N_{Se} = N_{total} \cdot \pi
\quad\text{and}\quad
N_{Sp} = N_{total} \cdot (1 - \pi)
$$

For sensitivity:

$$
z_{\beta_{Se}} = \frac{\sqrt{N_{Se}}|Se_1 - Se_0|  - z_{1-\alpha}\sqrt{Se_0(1-Se_0)}}{\sqrt{Se_1(1-Se_1)}}
$$

$$
Power_{Se} = \Phi(z_{\beta_{Se}})
$$

Similar for Specificity:

$$
z_{\beta_{Sp}} = \frac{\sqrt{N_{Sp}}|Sp_1 - Sp_0| - z_{1-\alpha}\sqrt{Sp_0(1-Sp_0)}}{\sqrt{Sp_1(1-Sp_1)}}
$$

$$
Power_{Sp} = \Phi(z_{\beta_{Sp}})
$$

The resulting overall power is

$$
Power_{overall} = Power_{Se} \times Power_{Sp}
$$

**Remark**  
Sample size determination and power evaluation rely on the same
asymptotic normal approximation framework, ensuring mathematical
consistency. With `exact = TRUE` the power is the exact binomial
rejection probability instead. It is a sawtooth function of $N$, so the
sample size is found by a scan over integer $N$ and `exact_rule` selects
either the smallest $N$ reaching the target (`"first"`) or the smallest
$N$ from which the power stays at or above the target (`"stays_above"`).

### Comparative designs

`power.copr.comp.test()` handles studies comparing an experimental test
with a comparator, for sensitivity and specificity, with superiority or
non-inferiority hypotheses (`type.se`, `type.sp`, margins `delta.se`,
`delta.sp`).

- `design = "unpaired"`: subjects are randomised to one of the two
  tests. The variance under the null is set by `variance` (`"null"`,
  `"pooled"` or `"unpooled"`).
- `design = "paired"`: all subjects receive both tests. The discordance
  probabilities `psi.d` and `psi.nd` must be given. Superiority uses the
  Miettinen formula and non-inferiority the Liu score method.

The significance level `alpha` is one-sided (default 0.025) as for the
single-test design. The overall power is again the product of the two
endpoint powers.

### Conventional sample size

`ss.copr.conv()` implements the conventional three-step approach of
Stark et al. (2022). It calculates the number of diseased subjects
needed for the sensitivity and the number of non-diseased subjects
needed for the specificity, each with its own endpoint power. These are
divided by the prevalence and by one minus the prevalence, and the
maximum is taken. The endpoint with the smaller requirement is
overpowered, so the overall power is typically above the target.

### Empirical (simulated) power

`power.copr.emp()` and `power.copr.comp.emp()` give vectorised Monte
Carlo estimates of the overall power for the single-test and the
unpaired comparative design. They are used to check the analytic results
and are not interchangeable, because the comparator is estimated in the
unpaired design.

## Example

### Required Sample Size

``` r

res <- power.copr.test(
  power = 0.80,
  prev = 0.5,
  se.e = 0.85, se.min = 0.70,
  sp.e = 0.85, sp.min = 0.65,
  alpha = 0.025
)

res
#> 
#>      Single-test design (sample size calculation)
#>               Wald normal approximation 
#> 
#>               n = 134
#>           n_dis = 67
#>        n_nondis = 67
#>            prev = 0.5
#>       sig.level = 0.025
#>     alternative = one.sided
#>           power = 0.8018057
#>        power_se = 0.8220385
#>        power_sp = 0.9753871
#>       ci_method = wald
#>           exact = FALSE
#>    target_power = 0.8
#> 
#> NOTE: n is the total sample size.
```

### Achieved Power for Fixed Sample Size

``` r

pow <- power.copr.test(
  N = res$n,
  prev = 0.5,
  se.e = 0.85, se.min = 0.70,
  sp.e = 0.85, sp.min = 0.65,
  alpha = 0.025
)
pow
#> 
#>      Single-test design (power calculation)
#>            Wald normal approximation 
#> 
#>               n = 134
#>           n_dis = 67
#>        n_nondis = 67
#>            prev = 0.5
#>       sig.level = 0.025
#>     alternative = one.sided
#>           power = 0.8018057
#>        power_se = 0.8220385
#>        power_sp = 0.9753871
#>       ci_method = wald
#>           exact = FALSE
#> 
#> NOTE: n is the total sample size.
```

### Comparative Design

``` r
power.copr.comp.test(
  power = 0.80,
  prev = 0.6,
  se.c = 0.75, se.e = 0.85,
  sp.c = 0.65, sp.e = 0.75
)
#> 
#>      Comparative unpaired design (sample size calculation)
#>                    Wald normal approximation 
#> 
#>               n = 909
#>           n_dis = 545.4
#>        n_nondis = 363.6
#>            prev = 0.6
#>       sig.level = 0.025
#>     alternative = one.sided
#>           power = 0.8003529
#>        power_se = 0.9784401
#>        power_sp = 0.8179886
#>    target_power = 0.8
#>         type.se = sup
#>         type.sp = sup
#> 
#> NOTE: n is the total sample size.
```

## References

- Stark, M., & Zapf, A. (2020). Sample size calculation and
  re-estimation based on the prevalence in a single-arm confirmatory
  diagnostic accuracy study. *Statistical Methods in Medical Research*,
  29(10), 2958–2971.

- Obuchowski, N. A. (1998). Sample size calculations in studies of test
  accuracy. *Statistical Methods in Medical Research*, 7(4), 371–392.

- U.S. Food and Drug Administration. (2022). *Multiple Endpoints in
  Clinical Trials: Guidance for Industry*. Silver Spring, MD: FDA.
  Available at: <https://www.fda.gov/media/162416/download>
