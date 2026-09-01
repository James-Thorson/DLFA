
# setwd( R'(C:\Users\jtuth\Documents\GitHub\DLFA\simulation)' )

library(dsem)
library(ggplot2)

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
  F1 -> area1, 0, area1_0, 0.707
  F1 -> area2, 0, area2_0, 0
  F1 -> area3, 0, area3_0, 0

  F1 -> area1, 1, area1_1, 0
  F1 -> area2, 1, area2_1, 0.707
  F1 -> area3, 1, area3_1, 0

  F1 -> area1, 2, area1_2, 0
  F1 -> area2, 2, area2_2, 0
  F1 -> area3, 2, area3_2, 0.707

  F1 <-> F1, 0, sd_F1, 1
  area1 <-> area1, 0, sd_area1, 0.707
  area2 <-> area2, 0, sd_area2, 0.707
  area3 <-> area3, 0, sd_area3, 0.707
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

sims = simulate(
  simulator,
  nsim = nr,
  variance = "random",
  resimulate_gmrf = TRUE
)

fit_r = vector( "list", length = nr )
for(ri in seq_len(nr)){
  fit_r[[ri]] = dsem(
    tsdata = sims[[ri]],
    sem = dynamics
  )
}

par_true = simulator$obj$par
true_df <- data.frame(
  par = factor(1:13),
  par_true = as.numeric(par_true)   # drop the names attribute
)

grid = expand.grid(
  par = seq_along(par_true),
  rep = seq_len(nr)
)
grid$parhat = as.vector(sapply( fit_r, \(x) x$opt$par ))

ggplot(grid, aes(x = factor(par), y = parhat)) +
  geom_violin(fill = "grey80", trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(
    data = true_df, aes(x = par, y = par_true),
    shape = "-", size = 10, color = "red"
  ) +
  labs(x = "Parameter", y = "Value")
ggsave( file = "simulation.png", width = 5, height = 4 )
