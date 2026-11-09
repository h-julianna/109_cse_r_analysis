# Bayesian Analysis and Precision Simulation for the Monetary Modulation of Conflict Adaptation

This repository contains the complete analysis pipeline, the parameter approximation and precision simulation used in the manuscript:

> **Probing Monetary Gain and Loss as Affective Signals for Adaptive Control**

The repository contains three main scripts:

| File | Purpose |
|------|---------|
| `109_cse_analysis_v2.Rmd` | Complete preprocessing, descriptive analyses, Bayesian hierarchical model fitting, model diagnostics, posterior summaries, and figures reported in the manuscript. |
 `109_precision_approximation.Rmd` | A markdown script that uses the pilot parameters to approximate the required N for the target HDI widths. |
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

Two models are fitted and compared with LOO:

- **`mu_only`**: predictors on `mu` only; `sigma` and `beta` (the exponential tail, tau) are intercept-only.
- **`mu_tau`**: the same predictors on `mu` and on `beta`; `beta` has a by-participant random intercept only.


### Contrast coding (`model_data`)

- `curr_cong`: incongruent = +0.5, congruent = −0.5
- `prev_cong`: incongruent = −0.5, congruent = +0.5
- `c_loss` / `c_gain`: non-orthogonal (2/3, −1/3, −1/3) sum-to-zero contrasts on the
  3-level monetary condition (`prev_color`). Despite non-orthogonality, because both
  contrasts are entered simultaneously, each three-way coefficient reduces to a clean
  loss-vs-neutral or gain-vs-neutral comparison in raw ms (see in-script comments for the
  derivation).

### Formulas

```r
# mu_only
exg_formula = bf(
  rt ~ curr_cong * prev_cong * c_loss + curr_cong * prev_cong * c_gain +
       (curr_cong * prev_cong * c_loss + curr_cong * prev_cong * c_gain | subj_code),
  family = exgaussian(link = "identity")
)

# mu_tau
exg_tau_formula = bf(
  rt ~ curr_cong * prev_cong * c_loss + curr_cong * prev_cong * c_gain +
       (curr_cong * prev_cong * c_loss + curr_cong * prev_cong * c_gain | subj_code),
  beta ~ curr_cong * prev_cong * c_loss + curr_cong * prev_cong * c_gain +
       (1 | subj_code),
  family = exgaussian(link = "identity", link_beta = "identity")
)
```

Note this is **not** a full three-way crossing of all predictors — `c_loss` and `c_gain`
each get their own three-way interaction with `curr_cong:prev_cong`, but `c_loss:c_gain`
and any four-way terms are excluded.

### Coefficient → hypothesis mapping

H1–H4 are tested on the `mu` coefficients. In brms, `mu` is already the effect on the
mean of the whole ExGaussian (mu = xi + beta), so `mu` and `beta` coefficients are **not**
summed. This holds for both models and the names are identical.

### Coefficient → hypothesis mapping

| Coefficient | Hypothesis | Interpretation |
|---|---|---|
| `b_curr_cong` | H1 (CE) | incongruent − congruent |
| `b_curr_cong:prev_cong` | H2 (CSE) | CSE interaction term |
| `b_curr_cong:prev_cong:c_loss` | H3 | CSE modulation, loss vs. neutral |
| `b_curr_cong:prev_cong:c_gain` | H4 | CSE modulation, gain vs. neutral |

In `mu_tau`, the corresponding `b_beta_*` coefficients describe effects on the ExGaussian
tail. They are not part of the H1–H4 inference.

### Priors

All on the raw-ms scale (identity link).

- Shared (both models): `Intercept ~ N(500,100)`; `curr_cong ~ N(0,50)`;
  `curr_cong:prev_cong ~ N(0,25)`; both three-way terms ~ `N(0,12.5)`; other fixed effects
  fall back to `N(0,50)`; `sd ~ N(0,50)`; `sigma ~ N(0,100)`.
- `mu_only`: `beta ~ N(0,100)` (intercept).
- `mu_tau`: `N(100,30)`; `N(0,20)`;
  `N(0,25)`.

### Sampling

4 chains, 6000 iterations (3000 warmup), `adapt_delta = 0.95` (raised from default to fix
low ESS/high Rhat on the intercept). The fits are cached as `fit_exg_model.rds` (`mu_only`)
and `fit_exg_tau_model.rds` (`mu_tau`), with LOO objects in `loo_exg.rds` and
`loo_exg_tau.rds`. These files are too large for GitHub and are hosted on Zenodo. Download them into the repository root to skip refitting; if a file is
missing, the Rmd refits that model.

## Model comparison

`mu_only` and `mu_tau` are compared with `loo::loo_compare()` (elpd difference and SE).
The better-fitting model is used for the confirmatory H1–H4 inference below. 

## Diagnostics

For both models: `summary()` (Rhat < 1.01, ESS > 400 target), trace plots for all H1–H4
coefficients, and a posterior predictive density overlay check. For `mu_tau`, also check
the `beta` coefficients.

## Inference: precision + ROPE

For each of H1–H4:

- **Precision (stopping rule):** 95% HDI width compared against a pre-specified target
  (H1 = 30ms, H2 = 20ms, H3 = 10ms, H4 = 10ms).
- **Inference (decision rule):** `bayestestR::equivalence_test()` against a ROPE of
  ±5ms, producing "effect present" / "practical null" / "inconclusive."

These two criteria are computed and reported separately, then combined into a single
`results` summary table (HDI bounds, width, target, precision met y/n, ROPE decision).
Each hypothesis also gets a `plotPost()` visualization with the HDI and ROPE marked.

## Manipulation check

We aggregated the three manipulation-check scores / trial type, to yield a single valence and arousal score per participant for each monetary condition. We compputed directional paired difference scores, and fitted an **intercept-only Gaussian model in brms.** Inference followed the same framework as the primary hypothesis. 

## Outputs

- `processed_pilot_data.csv` — cleaned analysis dataset
- `fit_exg_model.rds` — cached brms model fit, mu_only
- `fit_exg_tau_model.rds` — cached brms model fit, mu_tau
- Console/inline: descriptive stats, diagnostic plots, H1–H4 posterior summaries and
  precision/ROPE results table

# Precision simulation (`109_precision_simulation.R`)

Estimates the sample size needed for the 95% HDI width of H1–H4 to fall below the
targets (H1 = 30 ms, H2 = 20 ms, H3 = 8 ms, H4 = 8 ms), using the `mu_tau` ExGaussian
model. Targets are applied at the summary step, so they can be changed without refitting.

1. **True effects (mean RT):** CE = 50 ms, CSE = 25 ms, loss and gain modulation =
   12.5 ms each, with no effects on tau. These values come from the literature, **not**
   from the pilot. Only nuisance parameters (random-effect SDs, `sigma`, tau intercept
   and its by-participant SD) are taken from the pilot fit (`fit_exg_tau_model.rds`).
2. **Simulate:** fully crossed trial-level data (2 × 2 × 3 design, 37 trials per cell).
3. **Fit:** a simplified `mu_tau` model (full fixed effects on `mu` and `beta`; `mu`
   random effects for five terms only, uncorrelated; random intercept on `beta`),
   compiled once and reused via `update()`.
4. **Measure:** 95% HDI widths for H1–H4, bias, and sampler diagnostics.
5. **Run:** N = 500, 20 replicates.

# Software

R with `brms`, `Stan`, `bayestestR`, `posterior`, `dplyr`, `tidyr`, `purrr`, `MASS`.

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
It needs only the fitted pilot model `fit_exg_tau_model.rds`, which is hosted on Zenodo. Download it into the repository root.

The script reads the nuisance parameters from this fit (random-effect SDs, `sigma`,
and the tau intercept and its by-participant SD). The effect sizes used for data
generation (CE = 50 ms, CSE = 25 ms, loss and gain modulation = 12.5 ms) are set
explicitly in the script and are **not** estimated from the pilot.

Compiled Stan models are platform-specific, so the script recompiles the model on your
machine automatically (a short `update()` call with `recompile = TRUE`). This does not
refit the model, and no manual recompilation is needed.

Each replicate runs as one SLURM array task (`sbatch --array=1-20`); the script stops if
`SLURM_ARRAY_TASK_ID` is not set. To run a single replicate locally, set it first, e.g.
`Sys.setenv(SLURM_ARRAY_TASK_ID = 1)`.

# Notes

To reduce computation time, the simulation script

- reuses the compiled Stan model,
- fits shorter MCMC chains than the final analysis,
- evaluates multiple candidate sample sizes through repeated simulation.

---