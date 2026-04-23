# Analysis Code — README

## Overview

This repository contains the analysis code for the study, written in **R Markdown**.

---

## Data Preprocessing

Before statistical analyses, the following preprocessing steps are applied to the raw data:

- First trials of each block and all practice trials are excluded
- Incorrect trials and trials immediately preceded by an incorrect trial are excluded
- Trials with reaction times faster than **150 ms** are excluded
- Participants who did not fully complete the experiment are excluded
- Participants whose overall accuracy fell below **60%** are excluded

## Statistical Analyses

### Bayesian Tests
In addition to frequentist tests, Bayesian analyses are conducted using **Dienes's Bayes Factor calculator** (`old_Bf` function; Dienes, 2016). This function computes the relative likelihood of the data under the theoretical model and the null model via numerical integration, and returns their ratio as a Bayes Factor (BF).

#### Required inputs from the sampling distribution:
| Parameter | Description |
|---|---|
| `sd` | Standard error of the estimate |
| `obtained` | The observed effect (mean or mean difference) |
| `dfdata` | Degrees of freedom from the t-test |

#### Prior parameters for the theoretical model:
| Parameter | Description |
|---|---|
| `meanoftheory` | Mean of the theoretical prior distribution |
| `sdtheory` | Standard deviation of the theoretical prior |
| `dftheory` | Degrees of freedom of the theoretical distribution (set to `10^10` for a normal approximation) |
| `tail` | `1` for a one-tailed (half-normal) prior, `2` for two-tailed |

## Robustness Regions

To ensure that Bayes Factor conclusions are not driven by the specific prior width chosen, **robustness regions (RR)** are computed by systematically varying the `sdtheory` parameter across a range of values.

This is implemented using a custom loop that calls `old_Bf` iteratively:

**How this works:** `h1_range` defines a sequence of candidate prior standard deviations. For each value, `old_Bf` is called with that `sdtheory`, and the resulting BF is recorded. The output is a two-column matrix of `sdtheory` and `BF` values. The robustness region is then extracted as the range of `sdtheory` values for which BF > 3 (evidence in favour of the alternative), giving a lower and upper boundary. 

## Interpretation Thresholds

| BF | Interpretation |
|---|---|
| BF > 3 | Moderate evidence for the alternative hypothesis |
| 1/3 < BF < 3 | Inconclusive / insensitive |
| BF < 1/3 | Moderate evidence for the null hypothesis |




