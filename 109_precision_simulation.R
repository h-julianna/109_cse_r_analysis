# ------------------------------------------------------------------------------
# PRECISION-BY-SIMULATION for the hierarchical ExGaussian CE/CSE model
#
# The loop:
#   1. Hypothesize a "true" data-generating world     -> true_effects / noise
#   2. Simulate one dataset at candidate N            -> simulate_dataset()
#   3. Fit the SAME model used in the real analysis   -> goal_achieved_for_sample()
#   4. Tally whether each HDI width beat its target   -> goal_achieved_for_sample()
#   5. Repeat many times per N, sweep over N          -> main loop at bottom
#
# Requires: the pilot model object `fit_exg` already fit in 109_cse_analysis_v2.Rmd
# (this script assumes it is run in the same R session / project, after that
# model has been fit and cached).
#
# True_effects are literature-anchored (CE = 50ms, CSE = 25ms, loss/gain
# modulation = 12.5ms each).
# The pilot fit is still used for the nuisance parameters (random-effect
# SDs/correlation, sigma, beta), since those aren't the effects being
# powered on and the pilot is a reasonable source for typical noise levels.
# ------------------------------------------------------------------------------

library(tidyverse)
library(brms)
library(bayestestR)
library(MASS)   # for mvrnorm — note: MASS::select clashes with dplyr::select,
                # so we always call dplyr::select explicitly below if needed

# Faster compiler
library(cmdstanr)

# cmdstanr::check_cmdstan_toolchain()
# cmdstanr::install_cmdstan()
# cmdstanr::set_cmdstan_path("C:/Users/diano/.cmdstan/cmdstan-2.39.0")
# cmdstan_path()
 
# Parallelization and within-chain threading
library(parallel)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

parallel::detectCores(logical = FALSE)
parallel::detectCores(logical = TRUE)

nCores = detectCores()

model_data = read_csv("model_data.csv")

exg_priors = c(
  prior(normal(500, 100), class = Intercept),
  prior(normal(0, 50),   class = b, coef = "curr_cong"),
  prior(normal(0, 25),   class = b, coef = "curr_cong:prev_cong"),
  prior(normal(0, 12.5), class = b, coef = "curr_cong:prev_cong:c_loss"),
  prior(normal(0, 12.5), class = b, coef = "curr_cong:prev_cong:c_gain"),
  prior(normal(0, 50),   class = b),  # fallback for every other main effect
  prior(normal(0, 50),    class = sd),
  prior(normal(0, 100),   class = sigma),
  prior(normal(0, 100),   class = beta)
)

exg_formula = bf(rt ~ curr_cong * prev_cong * c_loss + curr_cong * prev_cong * c_gain +
                   (curr_cong | subj_code), family = exgaussian (link = "identity"))

true_effects = list(
  intercept = 500,   # plausible grand-mean RT in ms; matches the model's Intercept prior mean
  ce        = 50,    # H1: congruency effect
  cse       = 25,    # H2: congruency sequence effect
  loss_mod  = 12.5,  # H3: loss_CSE vs. neutral_CSE 
  gain_mod  = 12.5   # H4: gain_CSE vs. neutral_CSE 
)

# Precision targets.
target_widths = list(H1 = 30, H2 = 20, H3 = 10, H4 = 10)

n_trials_per_cell = 37   # needs to be defined before this test (same value used in the sweep below)

### Fitting the model


fit_exg = brm(
  formula = exg_formula,
  data    = model_data,
  prior   = exg_priors,
  chains  = 2,
  cores   = 2,
  iter    = 4000,
  warmup  = 2000,
  control = list(adapt_delta = 0.90), # up from the default 0.8: prior fit showed relatively low ESS/high Rhat on the intercept only; smaller steps improve mixing here.
  max_treedepth = 12, 
  threads = threading(4),
  seed    = 109,
  backend = "cmdstanr"
  #file    = "fit_exg_model_test"  # caches the fit; delete the file to refit
)
summary(fit_exg)
class(fit_exg$fit)
fit_exg$backend
prior_summary(fit_exg)

# ------------------------------------------------------------------------------
# 1. HYPOTHESIZED "TRUE" WORLD
# ------------------------------------------------------------------------------
# Effect sizes: literature-anchored, NOT pulled from the pilot posterior.


# Nuisance parameters (between-subject variability + trial-level noise): these
# aren't hypotheses being tested, so anchoring them on the pilot fit is fine 

# fit_exg = readRDS("fit_exg_model.rds")

# One-time forced recompile for the given machine. The .rds file bundles a
# compiled Stan binary that is platform-specific. It will not run on a
# different machine/OS/compiler than the one that originally fit it, even
# though the R object itself loads fine.
# recompile = TRUE rebuilds the binary for this platform; chains=1/iter=10
# keeps this fast since only a valid compiled
# program here is need, not a real fit.

# fit_exg = update(fit_exg, backend = "cmdstanr", recompile = TRUE, chains = 2, cores = 2, max_treedepth = 14, iter = 100, warmup = 50)

pilot_ranef_sd  = VarCorr(fit_exg)$subj_code$sd[, "Estimate"]        
pilot_ranef_cor = VarCorr(fit_exg)$subj_code$cor[2, "Estimate", 1] 
pilot_sigma_hat = summary(fit_exg)$spec_pars["sigma", "Estimate"]
pilot_beta_hat  = summary(fit_exg)$spec_pars["beta",  "Estimate"]



# ------------------------------------------------------------------------------
# 2. SIMULATE ONE DATASET AT A GIVEN SAMPLE SIZE
# ------------------------------------------------------------------------------
# Fully-crossed 2 (curr_cong) x 2 (prev_cong) x 3 (prev_color) design, using
# the same contrast coding as 109_cse_analysis_v2.Rmd. rexgaussian()
# draws directly from the ExGaussian likelihood, so no separate "noise on top
# of mu" step is needed — mu already is the location parameter of the ExGaussian.

simulate_dataset = function(n_subj, n_trials_per_cell, true_effects, ranef_sd, ranef_cor, sigma, beta)  {
  design = expand_grid(congruency  = c("congruent", "incongruent"), prev_congruency = c("congruent", "incongruent"), prev_color = c("loss", "gain", "neutral")
  ) %>%
    mutate(
      # Identical to the contrasts in 109_cse_analysis_v2.Rmd
      curr_cong = if_else(congruency      == "incongruent", 0.5, -0.5),
      prev_cong = if_else(prev_congruency == "incongruent", -0.5, 0.5),
      c_loss = case_when(
        prev_color == "loss"    ~  2/3,
        prev_color == "gain"    ~ -1/3,
        prev_color == "neutral" ~ -1/3
      ),
      c_gain = case_when(
        prev_color == "loss"    ~ -1/3,
        prev_color == "gain"    ~  2/3,
        prev_color == "neutral" ~ -1/3
      )
    )

  # Correlated by-subject random intercept & curr_cong slope, matching the
  # random-effect structure of `fit_exg` ( (curr_cong | subj_code) ).
  Sigma = matrix(c(ranef_sd[1]^2, ranef_cor * ranef_sd[1] * ranef_sd[2], ranef_cor * ranef_sd[1] * ranef_sd[2], ranef_sd[2]^2), nrow = 2)
  ranef_draws = MASS::mvrnorm(n_subj, mu = c(0, 0), Sigma = Sigma)
  colnames(ranef_draws) = c("u_intercept", "u_curr_cong")

  subj_df = tibble(subj_code = paste0("sim_", seq_len(n_subj))) %>%
    bind_cols(as_tibble(ranef_draws))

  subj_df %>%
    crossing(design) %>%
    slice(rep(seq_len(n()), each = n_trials_per_cell)) %>%
    mutate(
      mu = true_effects$intercept + u_intercept + (true_effects$ce + u_curr_cong) * curr_cong +
        true_effects$cse * curr_cong * prev_cong + true_effects$loss_mod * curr_cong * prev_cong * c_loss +
        true_effects$gain_mod * curr_cong * prev_cong * c_gain, 
      rt = brms::rexgaussian(n(), mu = mu, sigma = sigma, beta = beta),
      subj_code = as.factor(subj_code)
    )
}

# ------------------------------------------------------------------------------
# 3 + 4. FIT THE SIMULATED DATA AND TALLY GOAL ACHIEVEMENT
# ------------------------------------------------------------------------------
# Speed trick: `update()` reuses the ALREADY-COMPILED Stan model from fit_exg
# instead of recompiling from scratch for every simulated dataset. This is the
# single biggest time-saver available. iter/warmup are reduced from our final
# 6000/3000 for speed during the sweep — spot-check a handful of runs against
# full iterations before trusting the final numbers (bulk ESS still matters
# for a stable HDI-width estimate, just less critically than for the real fit).



goal_achieved_for_sample = function(sim_data, target_widths, iter = 700, warmup = 300) 
  {
    fit_sim = update(fit_exg, 
                     newdata = sim_data, 
                     iter = iter, 
                     warmup = warmup, 
                     chains = 1, 
                     cores = 1,
                     backend  = "cmdstanr",
                     adapt_delta = c(0.90),
                     max_treedepth = 15,
                     threads = threading(4),
                     seed = 109,
                     # stan_model_args = list(stanc_options = list("O1")),
                     refresh = 10
  )

  draws = as_draws_df(fit_sim)

  hdi_width = function(x) {
    h = hdi(x, ci = 0.95)
    h$CI_high - h$CI_low
  }

  bias = function(x, true_val) mean(x) - true_val

  w = list(
    H1 = hdi_width(draws$b_curr_cong),
    H2 = hdi_width(draws$`b_curr_cong:prev_cong`),
    H3 = hdi_width(draws$`b_curr_cong:prev_cong:c_loss`),
    H4 = hdi_width(draws$`b_curr_cong:prev_cong:c_gain`)
  )

  tibble(
    H1_width = w$H1, H1_hit = w$H1 < target_widths$H1,
    H2_width = w$H2, H2_hit = w$H2 < target_widths$H2,
    H3_width = w$H3, H3_hit = w$H3 < target_widths$H3,
    H4_width = w$H4, H4_hit = w$H4 < target_widths$H4,
    # Recovery diagnostic posterior mean vs. the true simulated value, in ms.
    H1_bias = bias(draws$b_curr_cong, true_effects$ce),
    H2_bias = bias(draws$`b_curr_cong:prev_cong`, true_effects$cse),
    H3_bias = bias(draws$`b_curr_cong:prev_cong:c_loss`, true_effects$loss_mod),
    H4_bias = bias(draws$`b_curr_cong:prev_cong:c_gain`, true_effects$gain_mod)
  )
}


# RUNTIME TEST 
# ------------------------------------------------------------------------------
# One-off timing check at the largest/slowest candidate N (150).


system.time({
  test_sim = simulate_dataset(
    n_subj = 50, n_trials_per_cell = n_trials_per_cell,
    true_effects = true_effects,
    ranef_sd = pilot_ranef_sd, ranef_cor = pilot_ranef_cor,
    sigma = pilot_sigma_hat, beta = pilot_beta_hat
    )
  })
 
 system.time({
  test_result = goal_achieved_for_sample(test_sim, target_widths)
 })
 
 summary(test_result$fit)

# ------------------------------------------------------------------------------
# 5. SWEEP OVER CANDIDATE SAMPLE SIZES
# ------------------------------------------------------------------------------

n_trials_per_cell = 37          # observed mean trials/cell/subject in processed_pilot_data.csv
n_sims_per_N      = 20        
candidate_Ns       = 200  

dir.create("simulation_results", showWarnings = FALSE)

power_results <- map_dfr(candidate_Ns, function(N) {
  
  message("Running N = ", N)
  
  for (i in seq_len(n_sims_per_N)) {
    
    outfile <- file.path(
      "simulation_results",
      sprintf("N_%03d_sim_%03d.rds", N, i)
    )
    
    # Skip completed simulations
    if (file.exists(outfile)) {
      message("  Skipping simulation ", i)
      next
    }
    
    message("  Simulation ", i)
    
    sim_data <- simulate_dataset(
      n_subj = N,
      n_trials_per_cell = n_trials_per_cell,
      true_effects = true_effects,
      ranef_sd = pilot_ranef_sd,
      ranef_cor = pilot_ranef_cor,
      sigma = pilot_sigma_hat,
      beta = pilot_beta_hat
    )
    
    result <- goal_achieved_for_sample(sim_data, target_widths)
    
    saveRDS(result, outfile)
  }
  
  # Read all completed simulations for this N
  sims <- map_dfr(
    list.files(
      "simulation_results",
      pattern = sprintf("^N_%03d_sim_.*\\.rds$", N),
      full.names = TRUE
    ),
    readRDS
  )
  
  sims %>%
    summarise(
      across(ends_with("_hit"), mean),
      across(ends_with("_bias"), mean)
    ) %>%
    mutate(N = N, .before = 1)
})



sim_results <- list.files("simulation_results/", full.names = TRUE) %>%
  sort() %>%
  map_dfr(readRDS) %>%
  mutate(
    N = rep(candidate_Ns, each = n_sims_per_N),
    .before = 1
  )

power_results <- sim_results %>%
  group_by(N) %>%
  summarise(
    across(ends_with("_hit"), mean, na.rm = TRUE),
    .groups = "drop"
  )

power_results %>%
  pivot_longer(
    cols = ends_with("_hit"),
    names_to = "hypothesis",
    values_to = "power"
  ) %>%
  mutate(
    hypothesis = stringr::str_remove(hypothesis, "_hit")
  ) %>%
  ggplot(aes(x = N, y = power, colour = hypothesis, group = hypothesis)) +
  geom_line(
    linewidth = 1,
    position = position_dodge(width = 3)
  ) +
  geom_point(
    size = 2,
    position = position_dodge(width = 3)
  ) +
  geom_hline(
    yintercept = 0.80,
    linetype = "dashed",
    colour = "grey40"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  labs(
    title = "Probability of achieving target HDI width by sample size",
    x = "Participants (N)",
    y = "Estimated power"
  ) +
  theme_minimal()

# ------------------------------------------------------------------------------
# Power curve plot
# ------------------------------------------------------------------------------


new_target_widths = list(H1 = 30, H2 = 20, H3 = 10, H4 = 10)

power_results_new_targets <- sim_results %>%
  mutate(
    H3_hit = H3_width < new_target_widths$H3,
    H4_hit = H4_width < new_target_widths$H4
  ) %>% 
  group_by(N) %>%
  summarise(
    across(ends_with("_hit"), mean, na.rm = TRUE),
    .groups = "drop"
  )


power_results_new_targets %>% 
pivot_longer(
  cols = ends_with("_hit"),
  names_to = "hypothesis",
  values_to = "power"
) %>%
  mutate(
    hypothesis = stringr::str_remove(hypothesis, "_hit")
  ) %>%
  ggplot(aes(x = N, y = power, colour = hypothesis, group = hypothesis)) +
  geom_line(
    linewidth = 1,
    position = position_dodge(width = 3)
  ) +
  geom_point(
    size = 2,
    position = position_dodge(width = 3)
  ) +
  geom_hline(
    yintercept = 0.80,
    linetype = "dashed",
    colour = "grey40"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  ) +
  labs(
    title = "Probability of achieving target HDI width by sample size",
    x = "Participants (N)",
    y = "Estimated power"
  ) +
  theme_minimal()
