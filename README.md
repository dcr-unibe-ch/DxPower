
<!-- README.md is generated from README.Rmd. Please edit that file -->

``` r
devtools::load_all()
#> ℹ Loading DxPower
#> Warning: le package 'testthat' a été compilé avec la version R 4.5.2
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
  - [`ss_opt_bin()`](#ss_opt_bin)
  - [`power_opt_bin()`](#power_opt_bin)
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

- One-sided superiority hypotheses

- Co-primary testing using an intersection–union framework

- Wald-type asymptotic normal approximations

Extensions to non-inferiority settings, continuous endpoints, and more
general co-primary objectives are planned.

## Statistical Framework

### Study Setting

We consider a confirmatory diagnostic accuracy study with:

- Sensitivity $Se$
- Specificity $Sp$
- Disease prevalence $\pi$
- Total sample size $N$

The number of diseased (cases) and non-diseased (control) subjects is:

$$N_{Se}=N\pi,  N_{Sp}=N(1-\pi)$$ Diseased and non-diseased subjects
form disjoint subgroups. Under standard binomial sampling assumptions,
the corresponding estimators of sensitivity and specificity are
asymptotically independent.

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

$$\hat{Se} \sim N\!\left( Se,\; \frac{Se(1-Se)}{N_{Se}} \right),\hat{Sp} \sim N\!\left( Sp,\; \frac{Sp(1-Sp)}{N_{Sp}} \right) $$

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

### `ss_opt_bin()`

Computes the required total sample size for a desired overall power and
true prevalence $\pi$ in a confirmatory diagnostic accuracy study with
co-primary endpoints.

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

### `power_opt_bin()`

Computes the achieved overall power for a fixed total sample size and
true prevalence $\pi$ in a confirmatory diagnostic accuracy study with
co-primary endpoints.

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
Both functions rely on the same asymptotic normal approximation
framework, ensuring mathematical consistency between sample size
determination and power evaluation.

## Example

### Required Sample Size

``` r

res <- ss_opt_bin(
  se1 = 0.85, sp1 = 0.85,
  se0 = 0.70, sp0 = 0.65,
  alpha = 0.025,
  prev = 0.5,
  pwr = 0.80
)

res
#> $pwr_se
#> [1] 0.820557
#> 
#> $pwr_sp
#> [1] 0.9749476
#> 
#> $N_se
#> [1] 66.77915
#> 
#> $N_sp
#> [1] 66.77915
#> 
#> $N_tot
#> [1] 133.5583
```

### Achieved Power for Fixed Sample Size

``` r

pow<-power_opt_bin(
  N = res$N_tot,
  se1 = 0.85,sp1 = 0.85, 
  se0 = 0.70,sp0 = 0.65,
  alpha = 0.025,
  prev = 0.5
)
pow
#> $overall
#> [1] 0.8
#> 
#> $pwr_se
#> [1] 0.820557
#> 
#> $pwr_sp
#> [1] 0.9749476
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
