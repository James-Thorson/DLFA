# Purpose: Sensitivity run to test sensitivity of fixed process sd, and how well stepwise selection recovers the lagged paths

library(dsem)
library(here)
library(dplyr)
library(parallel)

# ---- Fit the selected model ----------------------------------------------
# runs the data setup, saturated model, stepwise selection, and refit
source(here('R', 'DLFA.R'))

# ---- Value of fixed process SDs ------------------------------------------
sd_profile <- NULL
for (sd_fix in c(0.01, 0.1, 0.25, 0.5)) {
  fit_sd <- dsem(
    tsdata = tsdat,
    sem = gsub("NA, 0.01", paste0("NA, ", sd_fix), step_fit$model),
    family = fam,
    control = dsem_control(quiet = TRUE, newton_loops = 3)
  )
  s <- summary(fit_sd)
  s <- s[s$direction == 1 & s$lag > 0, c('path', 'lag', 'Estimate', 'Std_Error')]
  s$sd_fix <- sd_fix
  s$logLik <- -as.numeric(fit_sd$opt$objective)
  sd_profile <- rbind(sd_profile, s)
}
sd_profile

# ---- Parametric bootstrap simulation ------------------------------------------------
# simulate from the selected model (new factor, process error, and measurement
# error from the assessment SEs), set the factor to NA, refit the true model,
# and redo the stepwise selection (J = 7, T = 34, L = 3, same as the application)

# setup number of sims
n_boot <- 200
sims <- simulate(refit, nsim = n_boot, seed = 20261110, variance = 'random',
                 resimulate_gmrf = TRUE, fill_missing = TRUE)

# get true values
true_coefs <- s_step[s_step$direction == 1, c('path', 'lag', 'Estimate')]
names(true_coefs)[3] <- 'true'

# lagged paths that are in the selected model (the table from the last step of
# the selection has a 1 if a path is in the model, and a 0 if not)
last_step <- step_fit$step[[length(step_fit$step)]]
true_paths <- model_options[last_step['start', model_options] == 1]

# stepwise selection only needs the AIC, so skip the hessian steps
# (these give errors when a process SD goes to ~0)
ctrl_aic <- dsem_control(quiet = TRUE, getsd = FALSE, newton_loops = 0,
                         extra_convergence_checks = FALSE)

# setup bootstrap sim
run_boot <- function(b) {

  # make F1 a latent factor
  dat <- sims[[b]]
  F1_true <- as.numeric(dat[, 'F1'])
  dat[, 'F1'] <- NA

  # refit the true model. Same rule as the application: process SDs that
  # go to ~0 (< 0.05) are fixed at 0.01, and the model is fit again
  sem_b <- step_fit$model
  zero_sd <- NULL
  fit_b <- tryCatch(
    dsem(tsdata = dat, sem = sem_b, family = fam, control = ctrl_aic),
    error = function(e) NULL
  )
  if (!is.null(fit_b)) {
    sds <- summary(fit_b)
    sds <- sds[sds$direction == 2 & sds$parameter > 0, ] # estimated process SDs
    zero_sd <- sds$first[abs(sds$Estimate) < 0.05]
    for (r in zero_sd) sem_b <- paste0(sem_b, "\n", r, " <-> ", r, ", 0, NA, 0.01")
    fit_b <- tryCatch(
      dsem(tsdata = dat, sem = sem_b, family = fam,
           control = dsem_control(quiet = TRUE, newton_loops = 3)),
      error = function(e) NULL
    )
  }

  # redo the stepwise selection
  step_b <- tryCatch(
    stepwise_selection(model_options = model_options, model_shared = model_shared,
                       tsdata = dat, family = fam, quiet = TRUE, control = ctrl_aic),
    error = function(e) NULL
  )

  # extract out estimates
  est <- NULL
  if (!is.null(fit_b)) {

    # lag coefficients (sign of the factor is arbitrary, so flip to make these sum to > 0)
    coefs <- summary(fit_b)
    coefs <- coefs[coefs$direction == 1, ]
    flip <- sign(sum(coefs$Estimate))

    # estimated factor
    x_all <- fit_b$sdrep$par.random[names(fit_b$sdrep$par.random) == 'x_tj']
    F1_hat <- tail(x_all, Tn) * flip

    # put into df
    est <- data.frame(
      rep = b,
      path = coefs$path,
      lag = coefs$lag,
      estimate = coefs$Estimate * flip,
      se = coefs$Std_Error,
      pd_hess = isTRUE(fit_b$sdrep$pdHess),
      max_grad = max(abs(fit_b$obj$gr(fit_b$opt$par))),
      n_sd_fixed = length(zero_sd),
      factor_cor = cor(F1_hat, F1_true)
    )
  }

  # figure out what was kept (1 = lagged path is in the selected model, 0 = not)
  kept <- NULL
  if (!is.null(step_b)) {
    last_step <- step_b$step[[length(step_b$step)]]
    kept <- data.frame(
      rep = b,
      path = model_options,
      selected = as.numeric(last_step['start', model_options])
    )
  }

  list(est = est, kept = kept)
}

# run in parallel
boot <- mclapply(seq_len(n_boot), run_boot, mc.cores = detectCores() - 4)

# combine replicates
boot_est <- NULL
boot_kept <- NULL
for (b in seq_len(n_boot)) {
  boot_est <- rbind(boot_est, boot[[b]]$est)
  boot_kept <- rbind(boot_kept, boot[[b]]$kept)
}

# figure out retention of true and spurious lagged paths
length(unique(boot_kept$rep)) # number of replicates where selection finished
retention <- boot_kept %>%
  group_by(path) %>%
  summarize(retained = mean(selected)) %>%
  mutate(true_path = path %in% true_paths) %>%
  arrange(desc(true_path), desc(retained))
retention

# refits of the true model (converged = max gradient < 1e-3 and positive-definite hessian)
boot_est <- boot_est %>%
  left_join(true_coefs, by = c('path', 'lag')) %>%
  mutate(converged = pd_hess & max_grad < 1e-3)

# summary ... number converged, how often a process SD was fixed,
# and correlation between estimated and simulated factor
boot_rep <- boot_est %>% distinct(rep, converged, n_sd_fixed, factor_cor)
sum(boot_rep$converged)
mean(boot_rep$n_sd_fixed > 0)
median(boot_rep$factor_cor[boot_rep$converged])

# bias, RMSE, and coverage of 95% confidence intervals
boot_perf <- boot_est %>%
  filter(converged) %>%
  group_by(path, lag) %>%
  summarize(
    true = true[1],
    bias = mean(estimate - true),
    rmse = sqrt(mean((estimate - true)^2)),
    coverage = mean(abs(estimate - true) <= 1.96 * se),
    .groups = 'drop'
  )
boot_perf
