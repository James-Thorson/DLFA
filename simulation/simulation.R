
# setwd( R'(C:\Users\jtuth\Documents\GitHub\DLFA\simulation)' )

library(dsem)
library(ggplot2)
library(dplyr)

##############
# Parameterize simulator
##############

nt = 40

data = data.frame(
  area1 = rep(0, nt),
  area2 = rep(0, nt),
  area3 = rep(0, nt),
  F1 = rnorm(nt)
)

dynamics = "
  F1 -> area1, 0, area1_0, 1
  F1 -> area2, 0, area2_0, 0.5
  F1 -> area3, 0, area3_0, 0

  F1 -> area1, 1, area1_1, 0.5
  F1 -> area2, 1, area2_1, 1
  F1 -> area3, 1, area3_1, 0.5

  F1 -> area1, 2, area1_2, 0
  F1 -> area2, 2, area2_2, 0.5
  F1 -> area3, 2, area3_2, 1

  F1 <-> F1, 0, NA, 1
  area1 <-> area1, 0, sd_area1, 1
  area2 <-> area2, 0, sd_area2, 1
  area3 <-> area3, 0, sd_area3, 1
"

simulator = dsem(
  sem = dynamics,
  tsdata = ts(data),
  control = dsem_control(
    run_model = FALSE,
    quiet = TRUE
  )
)

##############
# Simulation experiment
##############

nr = 500

# seed has to be passed here (simulate.dsem resets it, and uses the clock if seed = NULL)
sims = simulate(
  simulator,
  nsim = nr,
  seed = 101,
  variance = "random",
  resimulate_gmrf = TRUE
)

par_true = simulator$obj$par
true_df <- data.frame(
  par = seq_along(par_true),
  par_true = as.numeric(par_true)
)
true_df$par_name = c(
  paste0("b", 1:3, "-0"),
  paste0("b", 1:3, "-1"),
  paste0("b", 1:3, "-2"),
  c( "sd1", "sd2", "sd3")
)

grid = NULL
converged = rep(FALSE, nr)
for(ri in seq_len(nr)){
  dat = sims[[ri]]
  dat[,'F1'] = NA # factor is latent (integrated out), not passed in as data
  fit = tryCatch(
    dsem(
      tsdata = dat,
      sem = dynamics,
      control = dsem_control(quiet = TRUE)
    ),
    error = function(e) NULL
  )
  if(is.null(fit)) next

  # converged = max gradient < 1e-3 and positive-definite hessian
  converged[ri] = isTRUE(fit$sdrep$pdHess) && max(abs(fit$obj$gr(fit$opt$par))) < 1e-3
  if(!converged[ri]) next

  # flip sign of factor so lag coefficients sum to > 0 (sign is arbitrary)
  p = fit$opt$par
  sgn = sign(sum(p[1:9]))
  p[1:9] = p[1:9] * sgn
  p[10:12] = abs(p[10:12])

  # estimated factor (conditional mode)
  re = fit$sdrep$par.random
  x_hat = tail(re[names(re) == 'x_tj'], nt) * sgn

  grid = rbind(grid, data.frame(
    rep = ri,
    par_name = true_df$par_name,
    par_true = true_df$par_true,
    parhat = p,
    se = sqrt(diag(fit$sdrep$cov.fixed)),
    factor_cor = cor(x_hat, sims[[ri]][,'F1'])
  ))
}
sum(converged)

# bias, RMSE, and coverage of 95% confidence intervals
perf = grid %>%
  group_by(par_name) %>%
  summarize(
    bias = mean(parhat - par_true),
    rmse = sqrt(mean((parhat - par_true)^2)),
    coverage = mean(abs(parhat - par_true) <= 1.96 * se, na.rm = TRUE)
  )
perf

# how often the largest coefficient is at the true dominant lag (area 1 = lag 0, area 2 = lag 1, area 3 = lag 2)
dom_lag = grid %>%
  filter(grepl("^b", par_name)) %>%
  mutate(area = as.numeric(substr(par_name, 2, 2)),
         lag = as.numeric(substr(par_name, 4, 4))) %>%
  group_by(rep, area) %>%
  summarize(lag_hat = lag[which.max(parhat)], .groups = 'drop') %>%
  group_by(area) %>%
  summarize(prop = mean(lag_hat == area - 1))
dom_lag

# correlation between estimated and simulated factor
median(grid$factor_cor)

# labels as beta_{area,lag} and sigma_area
labs_math = c(
  "b1-0" = "beta[1*','*0]", "b1-1" = "beta[1*','*1]", "b1-2" = "beta[1*','*2]",
  "b2-0" = "beta[2*','*0]", "b2-1" = "beta[2*','*1]", "b2-2" = "beta[2*','*2]",
  "b3-0" = "beta[3*','*0]", "b3-1" = "beta[3*','*1]", "b3-2" = "beta[3*','*2]",
  "sd1" = "sigma[1]", "sd2" = "sigma[2]", "sd3" = "sigma[3]"
)

ggplot(grid, aes(x = par_name, y = parhat)) +
  geom_hline(yintercept = 0, lty = 3, colour = "grey60") +
  geom_violin(fill = "grey80", trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(
    data = true_df, aes(x = par_name, y = par_true),
    shape = "-", size = 10, color = "red"
  ) +
  scale_x_discrete(labels = \(x) parse(text = labs_math[x])) +
  labs(x = "Parameter", y = "Estimate") +
  theme_bw()
ggsave( file = here::here("simulation", "simulation.png"), width = 7, height = 4 )
