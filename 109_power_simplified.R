# ==============================================================================
# PRECISION-BY-SIMULATION for the mu_tau ExGaussian model (study 109)
#
# Question: how many participants (N) are needed so that the 95% HDI width of
# the H1-H4 effects falls below the target widths (e.g. 8 ms for H3/H4)?
#
# Loop:
#   1. Define a "true" world (effects on MEAN RT, nuisance params from pilot fit)
#   2. Simulate one dataset at sample size N          -> simulate_dataset()
#   3. Fit the simplified mu_tau model              -> run_one()
#   4. Store HDI widths + diagnostics 
#   5. Repeat per N, sweep N; apply targets afterwards -> summarise_precision()
#
#
# Requires: fit_exg_tau_model.rds (original mu_tau pilot fit).
# Simplified mu random effects; original full fixed effects and beta model retained.
# ==============================================================================

library(dplyr)
library(brms)
library(bayestestR)
library(posterior)
library(purrr)
library(tidyr)


pilot_fit_file = "fit_exg_tau_model.rds"

out_dir = "simulation_results_mu_tau_simplified"

candidate_Ns = c(500)          # one sample size: 500 participants
n_sims_per_N = 20           # 20 independent SLURM array tasks
n_trials_per_cell = 37          # mean trials/cell/subject in the pilot

# Each SLURM array task handles exactly one replicate, numbered 1 to 20.
# Set SLURM_ARRAY_TASK_ID automatically with sbatch --array=1-20.
task_id = suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "")))
if (length(task_id) != 1L || is.na(task_id) ||
    !(task_id %in% seq_len(n_sims_per_N))) {
  stop("Run via sbatch --array=1-20 run_analysis.sh (SLURM_ARRAY_TASK_ID must be 1..20).")
}

# Skip completed tasks before loading the pilot or compiling the model.
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
N_to_run = candidate_Ns[[1L]]
outfile_to_run = file.path(out_dir,
  sprintf("N_%03d_sim_%03d.rds", N_to_run, task_id))
if (file.exists(outfile_to_run)) {
  message("Already finished: ", outfile_to_run)
  quit(save = "no", status = 0)
}
message("Starting simulation ", task_id, "/", n_sims_per_N,
        " at N=", N_to_run)

# Sampler settings for each simulated fit.
# Check n_div / max_rhat in the results before trusting.
sampler = list(chains = 2, cores = 2, iter = 2000, warmup = 1000,
                adapt_delta = 0.95, max_treedepth = 12)

# TRUE effects on MEAN RT (ms) 
true_total = c(H1 = 50, H2 = 25, H3 = 12.5, H4 = 12.5)

# TRUE effects on tau (ms). 0 = the whole effect sits on mu. The mean-RT effect
# stays at true_total either way (mu effect = total - tau effect).
true_tau = c(H1 = 0, H2 = 0, H3 = 0, H4 = 0)
mu_intercept = 500              # same as before; irrelevant for precision

# Precision targets (ms). Only used when summarising, so you can change them
# without re-running any fits.
target_widths = list(H1 = 30, H2 = 20, H3 = 10, H4 = 10)   # previous
new_target_widths = list(H1 = 30, H2 = 20, H3 = 8,  H4 = 8)    # new: 8 ms

hyp_terms = c(H1 = "curr_cong",
               H2 = "curr_cong:prev_cong",
               H3 = "curr_cong:prev_cong:c_loss",
               H4 = "curr_cong:prev_cong:c_gain")

# Same init as the real analysis (identity-link tau must start > 0)
init_fun = function() list(Intercept_beta = 100, b_beta = rep(0, 11))

# 1. parameters from the pilot mu_tau fit. 
# ------------------------------------------------------------------------------
fit_pilot = readRDS(pilot_fit_file)
pilot_draws = as_draws_df(fit_pilot)
pm = function(v) {
  stopifnot("parameter not found in pilot draws" = v %in% names(pilot_draws))
  mean(pilot_draws[[v]])
}

# Coefficient names, in the order model.matrix() produces them
design_formula = ~ curr_cong * prev_cong * c_loss + curr_cong * prev_cong * c_gain
coef_names = c("Intercept", "curr_cong", "prev_cong", "c_loss", "c_gain",
                "curr_cong:prev_cong", "curr_cong:c_loss", "prev_cong:c_loss",
                "curr_cong:c_gain", "prev_cong:c_gain",
                "curr_cong:prev_cong:c_loss", "curr_cong:prev_cong:c_gain")

# Reduced by-participant random-effects structure for mu.
# Keep the pilot posterior SDs for these five terms, but assume ZERO
# correlations, matching the || random-effects specification below.
random_terms = c("Intercept", "curr_cong", "curr_cong:prev_cong",
                 "curr_cong:prev_cong:c_loss", "curr_cong:prev_cong:c_gain")
stopifnot(all(random_terms %in% coef_names))
mu_sd = setNames(map_dbl(random_terms,
                         ~ pm(paste0("sd_subj_code__", .x))), random_terms)
Sigma_mu = diag(mu_sd^2, nrow = length(random_terms))
dimnames(Sigma_mu) = list(random_terms, random_terms)

# Coefficients (fixed effects) for mu and for beta (tau), same 12 names
beta_eff = setNames(rep(0, 12), coef_names)
beta_eff["Intercept"] = pm("b_beta_Intercept")
beta_eff[hyp_terms]   = true_tau

mu_eff = setNames(rep(0, 12), coef_names)
mu_eff["Intercept"] = mu_intercept
mu_eff[hyp_terms]   = true_total - true_tau

pars = list(
  mu_eff   = mu_eff,
  beta_eff = beta_eff,
  Sigma_mu = Sigma_mu,
  sd_tau_u = pm("sd_subj_code__beta_Intercept"),  # by-participant SD of tau
  sigma    = pm("sigma")
)
print(round(c(sigma = pars$sigma, tau_intercept = beta_eff[["Intercept"]],
              sd_tau_u = pars$sd_tau_u), 1))
print(round(sort(mu_sd, decreasing = TRUE), 1))   # selected by-participant SDs


# 2. simulate one dataset 
# ------------------------------------------------------------------------------
# Fully crossed 2 (curr_cong) x 2 (prev_cong) x 3 (prev_color) design

simulate_dataset = function(n_subj, n_trials_per_cell, pars, min_tau = 1) {
  design = expand_grid(congruency      = c("congruent", "incongruent"),
                        prev_congruency = c("congruent", "incongruent"),
                        prev_color      = c("loss", "gain", "neutral")) %>%
    mutate(
      curr_cong = if_else(congruency      == "incongruent",  0.5, -0.5),
      prev_cong = if_else(prev_congruency == "incongruent", -0.5,  0.5),
      c_loss = case_when(prev_color == "loss" ~ 2/3, TRUE ~ -1/3),
      c_gain = case_when(prev_color == "gain" ~ 2/3, TRUE ~ -1/3)
    )
  n_cells = nrow(design)

  X = model.matrix(design_formula, data = design)
  colnames(X) = sub("(Intercept)", "Intercept", colnames(X), fixed = TRUE)
  X = X[, coef_names]

  # Random effects exist only for the five selected terms.
  X_random = X[, random_terms, drop = FALSE]
  U = MASS::mvrnorm(n_subj, mu = rep(0, length(random_terms)),
                     Sigma = pars$Sigma_mu)
  U = matrix(U, nrow = n_subj, ncol = length(random_terms))
  u_tau = rnorm(n_subj, 0, pars$sd_tau_u)

  mu_mat   = sweep(U %*% t(X_random), 2, drop(X %*% pars$mu_eff), "+")      # subj x cell
  beta_mat = outer(u_tau, drop(X %*% pars$beta_eff), "+")            # subj x cell

  cells = design[rep(seq_len(n_cells), times = n_subj), ] %>%
    mutate(subj = rep(seq_len(n_subj), each = n_cells),
           mu   = as.vector(t(mu_mat)),
           beta = as.vector(t(beta_mat)))

  n_bad = sum(cells$beta < min_tau)
  if (n_bad > 0) warning(n_bad, " of ", nrow(cells), " subject-cells had tau < ", min_tau,
                         " ms; clipped.")

  cells[rep(seq_len(nrow(cells)), each = n_trials_per_cell), ] %>%
    mutate(rt        = brms::rexgaussian(n(), mu = mu, sigma = pars$sigma,
                                         beta = pmax(beta, min_tau)),
           subj_code = factor(sprintf("sim_%03d", subj))) %>%
    select(rt, curr_cong, prev_cong, c_loss, c_gain, subj_code)
}

# ------------------------------------------------------------------------------
# 3. COMPILE ONCE, THEN REFIT ON EACH SIMULATED DATASET
# ------------------------------------------------------------------------------
# We explicitly change the fitted formula: update(fit_pilot) by itself would
# silently keep the original, much larger random-effects model.
# Keep the original fixed effects for mu and beta (tau), and the original
# subject-level random intercept in beta.
new_model_formula = bf(
  rt ~ curr_cong * prev_cong * c_loss +
       curr_cong * prev_cong * c_gain +
       (1 + curr_cong + curr_cong:prev_cong +
        curr_cong:prev_cong:c_loss +
        curr_cong:prev_cong:c_gain || subj_code),
  beta ~ curr_cong * prev_cong * c_loss +
         curr_cong * prev_cong * c_gain + (1 | subj_code)
)

# Priors and family are inherited from fit_pilot. If the original pilot used
# coefficient-specific priors for random effects that are now absent, adjust
# those priors before fitting; brms will flag invalid prior specifications.
# The .rds can contain a platform-specific binary: compile once here.
set.seed(1)
fit_template = update(fit_pilot,
                      formula = new_model_formula,
                      newdata = simulate_dataset(20, 5, pars),
                      recompile = TRUE,
                      chains = 1, cores = 1, iter = 10, warmup = 5,
                      init = init_fun, refresh = 0)

hdi_width = function(x) { h = hdi(x, ci = 0.95); h$CI_high - h$CI_low }

run_one = function(N, sim_id, outfile = NULL) {
  sim_seed = 109000 + N * 100 + sim_id          # reproducible per (N, sim)
  set.seed(sim_seed)
  sim_data = simulate_dataset(N, n_trials_per_cell, pars)

  t0 = Sys.time()
  fit_sim = update(fit_template, newdata = sim_data, recompile = FALSE,
                    chains = sampler$chains, cores = sampler$cores,
                    iter = sampler$iter, warmup = sampler$warmup,
                    control = list(adapt_delta   = sampler$adapt_delta,
                                   max_treedepth = sampler$max_treedepth),
                    init = init_fun, seed = sim_seed, refresh = 0, silent = 2)
  runtime_min = as.numeric(difftime(Sys.time(), t0, units = "mins"))

  d  = as_draws_df(fit_sim)
  np = nuts_params(fit_sim)

  # mean-RT effect = mu coefficient + tau coefficient, per draw (option B)
  eff_total = imap(hyp_terms, ~ d[[paste0("b_", .x)]] + d[[paste0("b_beta_", .x)]])
  eff_mu    = imap(hyp_terms, ~ d[[paste0("b_", .x)]])

  key = c(paste0("b_", hyp_terms), paste0("b_beta_", hyp_terms),
           "b_Intercept", "b_beta_Intercept", "sigma")
  dg  = summarise_draws(subset_draws(d, variable = key), "rhat", "ess_bulk")

  res = tibble(
    N = N, sim = sim_id, runtime_min = runtime_min,
    max_rhat     = max(dg$rhat),
    min_ess_bulk = min(dg$ess_bulk),
    n_div        = sum(np$Value[np$Parameter == "divergent__"]),
    n_treedepth  = sum(np$Value[np$Parameter == "treedepth__"] >= sampler$max_treedepth)
  )
  for (h in names(hyp_terms)) {
    res[[paste0(h, "_width")]]    = hdi_width(eff_total[[h]])                 # primary
    res[[paste0(h, "_width_mu")]] = hdi_width(eff_mu[[h]])                    # mu coef only (reference)
    res[[paste0(h, "_bias")]]     = mean(eff_total[[h]]) - true_total[[h]]    # recovery check
  }
  if (!is.null(outfile)) saveRDS(res, outfile)
  res
}

# ------------------------------------------------------------------------------
# 4. SMOKE TEST (run once, interactively, before the sweep)
# ------------------------------------------------------------------------------
# Small N, just to confirm the pipeline runs end to end and to get a time estimate.
# Bias should be small relative to the width; n_div should be ~0.
#
# smoke = run_one(N = 40, sim_id = 0)
# print(smoke)

# ------------------------------------------------------------------------------
# 5. RUN ONLY THE SIMULATION ASSIGNED TO THIS ARRAY TASK
# ------------------------------------------------------------------------------
# A failure returns a non-zero exit code so SLURM records the task as failed.
# Results are uniquely named by N and simulation index.
run_one(N_to_run, task_id, outfile_to_run)
message("Finished: ", outfile_to_run)

# ------------------------------------------------------------------------------
# 6. SUMMARISE
# ------------------------------------------------------------------------------
load_results = function(dir = out_dir) {
  list.files(dir, pattern = "\\.rds$", full.names = TRUE) %>% map_dfr(readRDS)
}

summarise_precision = function(res, targets, assurance = 0.80) {
  long = res %>%
    select(N, sim, matches("^H[1-4]_width$")) %>%
    pivot_longer(matches("^H[1-4]_width$"), names_to = "hyp", values_to = "width") %>%
    mutate(hyp = sub("_width", "", hyp),
           target = unlist(targets)[hyp],
           hit = width < target)

  by_N = long %>%
    group_by(hyp, target, N) %>%
    summarise(n_sims = n(), median_width = median(width),
              hit_rate = mean(hit), .groups = "drop")

  # Approximation only: HDI width ~ 1/sqrt(N); validate on an N grid.
  # Using the `assurance` quantile of the widths gives the N at which about
  # `assurance` of simulated studies would reach the target. Works from a single N.
  n_needed = long %>%
    group_by(hyp, target, N) %>%
    summarise(w = quantile(width, assurance), .groups = "drop") %>%
    mutate(n_needed = N * (w / target)^2) %>%
    group_by(hyp, target) %>%
    summarise(n_needed_sqrt_rule = ceiling(median(n_needed)), .groups = "drop")

  # Joint assurance for the two most demanding hypotheses.
  joint_by_N = res %>%
    mutate(joint_hit_H3_H4 = H3_width < targets$H3 & H4_width < targets$H4) %>%
    group_by(N) %>%
    summarise(n_sims = n(), joint_hit_rate = mean(joint_hit_H3_H4),
              .groups = "drop")

  list(by_N = by_N, n_needed = n_needed, joint_by_N = joint_by_N)
}

plot_precision = function(by_N, title) {
  ggplot(by_N, aes(N, hit_rate, colour = hyp, group = hyp)) +
    geom_line(linewidth = 1, position = position_dodge(width = 3)) +
    geom_point(size = 2, position = position_dodge(width = 3)) +
    geom_hline(yintercept = 0.80, linetype = "dashed", colour = "grey40") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    labs(title = title, x = "Participants (N)",
         y = "P(95% HDI width < target)", colour = NULL) +
    theme_minimal()
}

# Aggregation is intentionally NOT executed inside each SLURM task.
# Once the 20 jobs finish, in an R session:
#   source("power_analysis_simplified.R")  # NOT recommended: triggers a simulation
# Instead use a separate short analysis script, for example:
#   library(dplyr); library(purrr); library(tidyr)
#   res = purrr::map_dfr(list.files("simulation_results_mu_tau_simplified",
#            pattern = "^N_500_sim_[0-9]+\\.rds$", full.names = TRUE), readRDS)
#   print(res)
# The summarise_precision() function above can be copied or sourced separately.
