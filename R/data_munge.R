library(here)
library(tidyverse)

# Read model outputs ----------------------------------------------------
sabie_rtmb_model <- readRDS(here("inputs", "5_rg_MltRel.RDS"))
nsei             <- readRDS(here("inputs", "actual_report.rds"))
nseidat          <- readRDS(here("inputs", "sablefish_data_19May2026.rds"))
nseisdrep        <- readRDS(here("inputs", "v23_3f_3s_2016_NEW_FINAL.RDS"))
wcsabie          <- readRDS(here("inputs", "wc_sabie.RDS"))

# Years to retain
years <- 1975:2025


# Helper to make a region x year data frame -----------------------------
make_ts <- function(x, regions, years, value_name) {

  stopifnot(
    nrow(x) == length(regions),
    ncol(x) == length(years)
  )

  expand_grid(
    Region = regions,
    Year = years
  ) %>%
    mutate(
      !!value_name := as.vector(t(x))
    )
}


# Region names ----------------------------------------------------------
regions_5 <- c("BS", "AI", "WGOA", "CGOA", "EGOA")

# RTMB model ------------------------------------------------------------
rtmb_year_idx <- which(1960:2025 %in% years)

rtmb_rec <- make_ts(
  sabie_rtmb_model$rep$Rec[, rtmb_year_idx, drop = FALSE],
  regions = regions_5,
  years = years,
  value_name = "Rec"
)

rtmb_ssb <- make_ts(
  sabie_rtmb_model$rep$SSB[, rtmb_year_idx, drop = FALSE],
  regions = regions_5,
  years = years,
  value_name = "SSB"
)

rtmb_recdev <- make_ts(
  sabie_rtmb_model$rep$ln_RecDevs[, rtmb_year_idx, drop = FALSE],
  regions = regions_5,
  years = years,
  value_name = "RecDevs"
)


# get rec devs and coerce
recdevs <- sqrt(diag(sabie_rtmb_model$sd_rep$cov.fixed)[names(diag(sabie_rtmb_model$sd_rep$cov.fixed)) == 'ln_RecDevs'])
recdevs_mat <- array(recdevs, dim = dim(sabie_rtmb_model$rep$ln_RecDevs))

rtmb_recdev_se <- make_ts(
  recdevs_mat[, rtmb_year_idx, drop = FALSE],
  regions = regions_5,
  years = years,
  value_name = "RecDevs_SE"
)



# SS3 model -------------------------------------------------------------
ss3_year_idx <- which(
  SPoRC::sgl_rg_wc_sablefish_data$years %in% years
)

ss3_years <- SPoRC::sgl_rg_wc_sablefish_data$years[ss3_year_idx]

ss3_rec <- rowSums(
  SPoRC::sgl_rg_wc_sablefish_data$ss3$NAA
)[ss3_year_idx]

ss3_ssb <- SPoRC::sgl_rg_wc_sablefish_data$ss3$SSB[
  ss3_year_idx
]

# Get recruitment deviations and SEs from the WC SPoRC model (age-0, indexed by year class)
sr  <- wcsabie$sdrep
map <- wcsabie$mapping$ln_RecDevs

# SEs for the 50 estimated devs, in map-level order
se_est <- sqrt(diag(sr$cov.fixed))[names(sr$par.fixed) == 'ln_RecDevs']

# scatter back into the full-length vector; fixed devs get NA
se_full <- rep(NA_real_, length(map))
se_full[!is.na(map)] <- se_est[as.integer(map[!is.na(map)])]

# year class y recruits at age-2 in year y + 2, so shift by 2 years to line up with
# the other series (all indexed by year of recruitment at age-2)
wc_idx <- match(ss3_years - 2, wcsabie$data$years)
ss3_recdev    <- as.vector(wcsabie$rep$ln_RecDevs)[wc_idx]
ss3_recdev_se <- se_full[wc_idx]

# SS3 is Alaska-wide, so give it a single region ------------------------
ss3_ts <- tibble(
  Region   = "WC",
  Year     = ss3_years,
  Rec      = ss3_rec,
  SSB      = ss3_ssb,
  RecDevs  = ss3_recdev,
  RecDevs_SE = ss3_recdev_se
)


# NSEI model ------------------------------------------------------------

# get recdevs and coerce
nsei_rec <- tibble(
  Region   = "CS",
  Year     = years,
  Rec      = nsei$pred_rec,
  SSB      = rowSums(nsei$spawn_biom)[-52],
  RecDevs  = nseisdrep$rep$par.fixed[names(nseisdrep$rep$par.fixed) == 'log_rec_devs'], # log rec devs (age-2), same as SEs
  RecDevs_SE =  sqrt(diag(nseisdrep$rep$cov.fixed)[names(diag(nseisdrep$rep$cov.fixed)) == 'log_rec_devs'])
)


# Combine ---------------------------------------------------------------
ts_df <- bind_rows(
  rtmb_rec %>%
    left_join(rtmb_ssb, by = c("Region", "Year")) %>%
    left_join(rtmb_recdev, by = c("Region", "Year")) %>%
    left_join(rtmb_recdev_se, by = c("Region", "Year")),
  ss3_ts,

  nsei_rec
) %>%
  mutate(
    Region = factor(
      Region,
      levels = c("BS", "AI", "WGOA", "CGOA", "EGOA", "WC", "CS")
    )
  ) %>%
  arrange(Region, Year)


# Write out --------------------------------------------------------------
write.csv(
  ts_df,
  here("outputs", "ts_df.csv"),
  row.names = FALSE
)
