library(dsem)
library(here)
library(dplyr)
library(ggplot2)
library(cowplot)

# ---- Data ----------------------------------------------------------------

long <- read.csv(here('outputs', 'ts_df.csv'))
long <- subset(long, Year >= 1990 & Year <= 2023)

wide <- tapply(
  long$RecDevs,
  INDEX = list(Year = long$Year, Region = long$Region),
  FUN = mean
)

# matching SE matrix, same shape / ordering as `wide`
wide_se <- tapply(
  long$RecDevs_SE,
  INDEX = list(Year = long$Year, Region = long$Region),
  FUN = mean
)

# standardize
wide_raw <- wide
wide <- scale(wide)
ctr <- attr(wide, "scaled:center")
sdv <- attr(wide, "scaled:scale")

regions <- colnames(wide)
yrs <- as.numeric(rownames(wide))
Tn <- length(yrs)
J <- length(regions)

reg_levels <- c('BS', 'AI', 'WGOA', 'CGOA', 'EGOA', 'CS', 'WC')
reg_levels <- reg_levels[reg_levels %in% regions]

tsdat <- ts(cbind(wide, F1 = NA), start = min(yrs))

# SD matrix with a column for the latent factor (NA there)
sd_mat <- cbind(wide_se, F1 = NA)

# create family
fam <- c(
  lapply(regions, function(r) dsem::gaussian_fixed_sd(sd = wide_se[, r], link = 'identity')),
  list(dsem::fixed())
)

names(fam) <- c(regions, "F1")

# ---- Saturated model -----------------------------------------------------
sem <- "
  F1 -> AI, 0, bAI_0
  F1 -> BS, 0, bBS_0
  F1 -> CGOA, 0, bCGOA_0
  F1 -> EGOA, 0, bEGOA_0
  F1 -> WGOA, 0, bWGOA_0
  F1 -> CS, 0, bCS_0
  F1 -> WC, 0, bWC_0

  F1 -> AI, 1, bAI_1
  F1 -> BS, 1, bBS_1
  F1 -> CGOA, 1, bCGOA_1
  F1 -> EGOA, 1, bEGOA_1
  F1 -> WGOA, 1, bWGOA_1
  F1 -> CS, 1, bCS_1
  F1 -> WC, 1, bWC_1

  F1 -> AI, 2, bAI_2
  F1 -> BS, 2, bBS_2
  F1 -> CGOA, 2, bCGOA_2
  F1 -> EGOA, 2, bEGOA_2
  F1 -> WGOA, 2, bWGOA_2
  F1 -> CS, 2, bCS_2
  F1 -> WC, 2, bWC_2

  F1 -> AI, 3, bAI_3
  F1 -> BS, 3, bBS_3
  F1 -> CGOA, 3, bCGOA_3
  F1 -> EGOA, 3, bEGOA_3
  F1 -> WGOA, 3, bWGOA_3
  F1 -> CS, 3, bCS_3
  F1 -> WC, 3, bWC_3

  F1 <-> F1, 0, NA, 1
"

# RUN WITH `getsd = FALSE` to see constraints needed
fit <- dsem(
  tsdata = tsdat,
  sem = sem,
  family = fam,
  control = dsem_control(
    trace = 1,
    getsd = FALSE
  )
)

summary(fit)

# adding constraints on process error
sem_constrained = paste0( sem,
"
  AI <-> AI, 0, NA, 0.01
  BS <-> BS, 0, NA, 0.01
")

# RUN again
fit <- dsem(
  tsdata = tsdat,
  sem = sem_constrained,
  family = fam,
  control = dsem_control(
    trace = 1,
    getsd = TRUE
  )
)


s_sat <- summary(fit)
s_sat

fit$sdrep

# ---- Stepwise path selection ---------------------------------------------
model_options <- c(
  "F1 -> AI, 1, bAI_1",
  "F1 -> BS, 1, bBS_1",
  "F1 -> CGOA, 1, bCGOA_1",
  "F1 -> EGOA, 1, bEGOA_1",
  "F1 -> WGOA, 1, bWGOA_1",
  "F1 -> CS, 1, bCS_1",
  "F1 -> WC, 1, bWC_1",
  "F1 -> AI, 2, bAI_2",
  "F1 -> BS, 2, bBS_2",
  "F1 -> CGOA, 2, bCGOA_2",
  "F1 -> EGOA, 2, bEGOA_2",
  "F1 -> WGOA, 2, bWGOA_2",
  "F1 -> CS, 2, bCS_2",
  "F1 -> WC, 2, bWC_2",
  "F1 -> AI, 3, bAI_3",
  "F1 -> BS, 3, bBS_3",
  "F1 -> CGOA, 3, bCGOA_3",
  "F1 -> EGOA, 3, bEGOA_3",
  "F1 -> WGOA, 3, bWGOA_3",
  "F1 -> CS, 3, bCS_3",
  "F1 -> WC, 3, bWC_3"
)

model_shared <- "
  F1 -> AI, 0, bAI_0
  F1 -> BS, 0, bBS_0
  F1 -> CGOA, 0, bCGOA_0
  F1 -> EGOA, 0, bEGOA_0
  F1 -> WGOA, 0, bWGOA_0
  F1 -> CS, 0, bCS_0
  F1 -> WC, 0, bWC_0

  AI <-> AI, 0, NA, 0.01
  BS <-> BS, 0, NA, 0.01

  F1 <-> F1, 0, NA, 1
"

# do step fit
step_fit <- dsem::stepwise_selection(
  model_options = model_options,
  model_shared  = model_shared,
  tsdata = tsdat,
  family = fam,
  quiet = TRUE
)

cat(step_fit$model)

# ---- Refit the selected model --------------------------------------------

refit <- dsem(
  tsdata = tsdat,
  sem = step_fit$model,
  family = fam,
  control = dsem_control(trace = 1, newton_loops = 3)
)

s_step <- summary(refit)
s_step

refit$sdrep
AIC(fit)
AIC(refit)

# ---- Extract the latent factor -------------------------------------------
x_all <- refit$sdrep$par.random[names(refit$sdrep$par.random) == 'x_tj']
v_all <- refit$sdrep$diag.cov.random[names(refit$sdrep$par.random) == 'x_tj']

x <- tail(x_all, Tn)

dlfa_trend <- data.frame(
  Year    = yrs,
  dlfa    = x,
  dlfa_se = sqrt(tail(v_all, Tn))
)

sd(dlfa_trend$dlfa)

# ---- Build the convolution -----------------------------------------------
conv <- sapply(regions, function(r) {
  br <- s_step[s_step$direction == 1 & s_step$second == r, c("lag", "Estimate")]
  if (nrow(br) == 0) return(rep(0, Tn))
  br <- br[order(br$lag), ]
  sapply(seq_along(x), function(t) {
    idx <- t - br$lag
    use <- idx >= 1                    # drop lags reaching before year 1
    sum(br$Estimate[use] * x[idx[use]])
  })
})
colnames(conv) <- regions

# ---- Plotting frame ------------------------------------------------------

plot_df <- data.frame(
  Year   = rep(yrs, times = J),
  Region = rep(regions, each = Tn),
  conv   = as.vector(conv),
  obs    = as.vector(wide)
) %>%
  left_join(dlfa_trend[, c('Year', 'dlfa')], by = 'Year') %>%
  mutate(Region = factor(Region, levels = reg_levels))

# ---- Convolution vs observed ---------------------------------------------

p_conv <- ggplot(plot_df, aes(x = Year)) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_line(aes(y = obs), lty = 2, lwd = 1.1, col = 'blue') +
  geom_line(aes(y = conv), lwd = 1.2) +
  facet_wrap(~Region) +
  labs(y = "Standardized recruitment deviation") +
  theme_bw(base_size = 15)

p_conv

# ---- Observed vs factor vs convolution -----------------------------------

p_comp <- ggplot(plot_df, aes(x = Year)) +
  geom_hline(yintercept = 0, linetype = "dotted",
             linewidth = 0.5, colour = "grey60") +
  geom_line(aes(y = obs,  colour = "Observed",    linetype = "Observed"),
            linewidth = 1.0) +
  geom_line(aes(y = dlfa, colour = "DLFA",        linetype = "DLFA"),
            linewidth = 1.1) +
  geom_line(aes(y = conv, colour = "Convolution", linetype = "Convolution"),
            linewidth = 1.1) +
  facet_wrap(~Region, ncol = 1) +
  scale_colour_manual(values = c(
    "Observed"    = "#0072B2",
    "DLFA"        = "#D55E00",
    "Convolution" = "#333333"
  )) +
  scale_linetype_manual(values = c(
    "Observed"    = "dashed",
    "DLFA"        = "solid",
    "Convolution" = "solid"
  )) +
  labs(x = 'Year', y = "Standardized recruitment deviation",
       colour = 'Type', linetype = 'Type') +
  theme_bw(base_size = 17) +
  theme(legend.position = 'top')

p_comp

# ---- Lag profiles: saturated vs selected ---------------------------------

get_coefs <- function(smry, label) {
  d <- smry[smry$direction == 1, c("second", "lag", "Estimate", "Std_Error")]
  names(d)[1] <- "Region"
  d$Model <- label
  d
}

coef_cmp <- rbind(
  get_coefs(s_sat,  "Saturated"),
  get_coefs(s_step, "Selected")
) %>%
  mutate(
    Region = factor(Region, levels = reg_levels),
    Model  = factor(Model, levels = c("Saturated", "Selected"))
  )

# fill dropped paths with exact zeros so gaps are visible
grid <- expand.grid(
  Region = factor(reg_levels, levels = reg_levels),
  lag    = 0:3,
  Model  = factor(c("Saturated", "Selected"),
                  levels = c("Saturated", "Selected"))
)

coef_cmp <- grid %>%
  left_join(coef_cmp, by = c("Region", "lag", "Model")) %>%
  mutate(dropped  = is.na(Estimate),
         Estimate = ifelse(dropped, 0, Estimate))

p_coef_cmp <- ggplot(coef_cmp %>% filter(Estimate != 0),
                     aes(x = lag, y = Estimate, colour = Model, shape = Model)) +
  geom_hline(yintercept = 0, lty = 3, colour = "grey60") +
  geom_pointrange(
    aes(ymin = Estimate - 1.96 * Std_Error,
        ymax = Estimate + 1.96 * Std_Error),
    position = position_dodge(width = 0.4)
  ) +
  facet_wrap(~Region, ncol = 1) +
  scale_colour_manual(values = c("Saturated" = "#999999",
                                 "Selected"  = "#000000")) +
  labs(x = "Lag (years)", y = "Coefficient") +
  theme_bw(base_size = 17) +
  theme(legend.position = "top")

p_coef_cmp

# ---- Combined figure -----------------------------------------------------

ggsave(
  filename = 'ts_lag_comp.png',
  plot = cowplot::plot_grid(p_comp, p_coef_cmp, ncol = 2,
                            label_size = 30,
                            labels = c("A", "B")),
  width = 14, height = 16, dpi = 300,
  path = here('figs')
)

# ---- Correlations: does the lag structure help? --------------------------
cor_check <- data.frame(
  Region     = regions,
  cor_factor = sapply(regions, function(r) cor(x, wide[, r], use = "complete.obs")),
  cor_conv   = sapply(regions, function(r) cor(conv[, r], wide[, r], use = "complete.obs"))
)
cor_check

