# Bayesian Analysis and Precision Simulation for the Monetary Modulation of Conflict Adaptation

This repository contains the complete analysis pipeline and parameter precision simulation used in the manuscript:

> **Probing Monetary Gain and Loss as Affective Signals for Adaptive Control**

The repository contains two main scripts:

| File | Purpose |
|------|---------|
| `109_cse_analysis_v2.Rmd` | Complete preprocessing, descriptive analyses, Bayesian hierarchical model fitting, model diagnostics, posterior summaries, and figures reported in the manuscript. |
| `109_precision_simulation.R` | Precision-by-simulation procedure used to determine the target sample size before final data collection. |

---

# Project overview

The study investigates whether monetary incentives modulate the **Congruency Sequence Effect (CSE)** in a prime-probe conflict task.

The project adopts a **Bayesian parameter estimation** framework. Data collection is guided by **posterior precision**, stopping once the target Highest Density Interval (HDI) widths for the focal parameters become sufficiently narrow.

The analysis is performed using a hierarchical Bayesian ExGaussian regression model implemented in **brms** (Stan backend).

---

# Repository structure

```
.
├── 109_cse_analysis_v2.Rmd
├── 109_precision_simulation.R
├── fit_exg_models.rds
├── raw_pilot_data.csv
└── README.md
```

---

# 109_cse_analysis_v2.Rmd

Analysis script for the prime-probe CSE study: cleans pilot data, fits a hierarchical
Bayesian ExGaussian model of reaction time, and evaluates H1–H4 using an HDI/ROPE
precision-estimation framework. H5–H8 (manipulation checks) are stubbed but not yet
implemented.

## Requirements

R 4.5.2, with packages: `tidyverse`, `brms`, `bayesplot`, `tidybayes`, `bayestestR`,
`posterior`, `see`, `bayesboot`.

## Input

`raw_pilot_data.csv` — raw trial-level export from the experiment software, expected in
the working directory.

## Pipeline

1. **Demographics** — extracts age/gender from embedded survey JSON responses.
2. **Practice trial flagging** — marks trials before the "Blokk 1 kezdődik" marker as
   practice and excludes them downstream.
3. **Previous-trial flagging** — computes `prev_congruency`, `prev_correct`, and
   `first_trial` via a 4-trial lag (accounts for trial structure where diagnostic trials
   are interleaved with inducer/monetary trials).
4. **Diagnostic trial flagging** — restricts to `task == "probe"` experimental trials;
   labels trials as `inducer` (monetary) or `diagnostic` (CE/CSE trials of interest).
5. **Monetary-monetary sequence exclusion** — drops trial pairs where a monetary trial
   directly follows another monetary trial (a randomization artifact that shouldn't occur
   by design).
6. **Accuracy filter** — excludes subjects with <60% overall accuracy; also drops
   `first_trial` rows (no valid previous trial).
7. **Correct-trials filter** — keeps only trials where both the current and previous
   trial were answered correctly, and RT ≥ 150ms.
8. **Column pruning** — drops housekeeping/unused columns, producing
   `processed_pilot_data.csv` (the main analysis dataset, written to disk).
9. **Previous-color flagging** — maps the previous trial's stimulus color to
   `loss`/`gain`/`neutral`, keeping only diagnostic trials with a defined previous color.

## Descriptives & plots

Includes an `dstats()` helper (min/max/mean/median/sd/skewness/kurtosis), RT distribution
checks (proportion of RTs > 1000ms), participant-level CE and CSE calculations,
and associated density/violin/interaction plots. These are exploratory and precede the
confirmatory model.

## Hierarchical ExGaussian model

### Contrast coding (`model_data`)

- `curr_cong`: incongruent = +0.5, congruent = −0.5
- `prev_cong`: incongruent = −0.5, congruent = +0.5
- `c_loss` / `c_gain`: non-orthogonal (2/3, −1/3, −1/3) sum-to-zero contrasts on the
  3-level monetary condition (`prev_color`). Despite non-orthogonality, because both
  contrasts are entered simultaneously, each three-way coefficient reduces to a clean
  loss-vs-neutral or gain-vs-neutral comparison in raw ms (see in-script comments for the
  derivation).

### Formula

```r
rt ~ curr_cong * prev_cong * c_loss + curr_cong * prev_cong * c_gain + (curr_cong | subj_code)
```

Family: `exgaussian(link = "identity")`. Note this is **not** a full three-way crossing of
all predictors — `c_loss` and `c_gain` each get their own three-way interaction with
`curr_cong:prev_cong`, but `c_loss:c_gain` and any four-way terms are excluded.

### Coefficient → hypothesis mapping

| Coefficient | Hypothesis | Interpretation |
|---|---|---|
| `b_curr_cong` | H1 (CE) | incongruent − congruent |
| `b_curr_cong:prev_cong` | H2 (CSE) | CSE interaction term |
| `b_curr_cong:prev_cong:c_loss` | H3 | CSE modulation, loss vs. neutral |
| `b_curr_cong:prev_cong:c_gain` | H4 | CSE modulation, gain vs. neutral |

### Priors

All on the raw-ms scale (identity link): `Intercept ~ N(500,100)`; `curr_cong ~ N(0,50)`;
`curr_cong:prev_cong ~ N(0,25)`; both three-way terms ~ `N(0,12.5)`; other fixed effects
fall back to `N(0,50)`; `sd ~ N(0,50)`; `sigma ~ N(0,100)`; `beta ~ N(0,100)`.

### Sampling

4 chains, 6000 iterations (3000 warmup), `adapt_delta = 0.95` (raised from default to fix
low ESS/high Rhat on the intercept). Cached to `fit_exg_model.rds` — delete this file to
force a refit.

## Diagnostics

`summary(fit_exg)` (Rhat < 1.01, ESS > 400 target), trace plots for all H1–H4 coefficients,
and a posterior predictive density overlay check.

## Inference: precision + ROPE

For each of H1–H4:

- **Precision (stopping rule):** 95% HDI width compared against a pre-specified target
  (H1 = 30ms, H2 = 20ms, H3 = 15ms, H4 = 15ms).
- **Inference (decision rule):** `bayestestR::equivalence_test()` against a ROPE of
  ±5ms, producing "effect present" / "practical null" / "inconclusive."

These two criteria are computed and reported separately, then combined into a single
`results` summary table (HDI bounds, width, target, precision met y/n, ROPE decision).
Each hypothesis also gets a `plotPost()` visualization with the HDI and ROPE marked.

## Manipulation check

We aggregated the three manipulation-check scores / trial type, to yield a single valence and arousal score per participant for each monetary condition. We compputed directional paired difference scores, and fitted an **intercept-only Gaussian model in brms.** Inference followed the same framework as the primary hypothesis. 

## Outputs

- `processed_pilot_data.csv` — cleaned analysis dataset
- `fit_exg_model.rds` — cached brms model fit
- Console/inline: descriptive stats, diagnostic plots, H1–H4 posterior summaries and
  precision/ROPE results table

# Precision simulation (`109_precision_simulation.R`)

## Step 1 — Specify a hypothetical "true" world

Expected effect sizes are anchored to previous literature while remaining intentionally conservative.

| Effect | Simulated value |
|---------|----------------|
| CE | 50 ms |
| CSE | 25 ms |
| Loss modulation | 12.5 ms |
| Gain modulation | 12.5 ms |

These values are **not estimated from the pilot data**.

Instead, the pilot model contributes only nuisance parameters such as

- participant variability,
- random-effect correlations,
- residual variance,
- ExGaussian parameters.

---

## Step 2 — Simulate datasets

For each candidate sample size the script generates fully crossed trial-level datasets matching the experimental design.

---

## Step 3 — Fit the analysis model

Every simulated dataset is analysed using **the exact same hierarchical model** employed in the final manuscript.

To improve computational efficiency,

`update()` is used to reuse the already-compiled Stan model instead of recompiling it for every simulation.

---

## Step 4 — Measure posterior precision

For each fitted model, the script computes the 95% HDI width for the four focal regression coefficients.

---

## Step 5 — Sweep across candidate sample sizes

The entire simulation is repeated many times for each candidate sample size.

The final output reports

- the proportion of successful simulations,
- average parameter bias,
- estimated probability of achieving the desired precision.

The script additionally generates a precision curve illustrating how the probability of meeting each HDI target changes as sample size increases.

---

# Software

The analyses were conducted in R using, among others,

- `brms`
- `Stan`
- `tidyverse`
- `bayestestR`
- `MASS`

---

# Reproducing the analyses

## Analysis

Run

```
109_cse_analysis_v2.Rmd
```

to

- fit the Bayesian hierarchical model,
- generate descriptive statistics,
- produce posterior summaries,
- create manuscript figures.

---

## Precision simulation

The precision simulation can be reproduced **without re-running the main analysis**.

The repository includes the fitted Bayesian model object from the pilot

```r
fit_exg_model.rds
```

which provides the nuisance parameter estimates (random-effect variances, correlations, and ExGaussian parameters) required for the simulations.

Because compiled Stan executables are platform-specific, the model should first be loaded and recompiled on the local machine:

```r
fit_exg = readRDS("fit_exg_model.rds")

fit_exg = update(
  fit_exg,
  recompile = TRUE,
  chains = 1,
  iter = 10,
  warmup = 5
)
```

The short update call does **not** refit the model for inference. It simply recompiles the underlying Stan program for the user's operating system, allowing the simulation script to reuse the compiled model efficiently.

The script uses the fitted model object only to obtain nuisance parameter estimates and to reuse the compiled Stan model during repeated simulation. The hypothesized effect sizes used for data generation (CE = 50 ms, CSE = 25 ms, loss modulation = 12.5 ms, gain modulation = 12.5 ms) are specified explicitly within the simulation script and are **not** estimated from the pilot model.

# Notes

To reduce computation time, the simulation script

- reuses the compiled Stan model,
- fits shorter MCMC chains than the final analysis,
- evaluates multiple candidate sample sizes through repeated simulation.

---