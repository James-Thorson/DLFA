library(here)
library(tidyverse)

# Read model outputs ----------------------------------------------------
sabie_rtmb_model <- readRDS(here("inputs", "5_rg_MltRel.RDS"))
nsei             <- readRDS(here("inputs", "actual_report.rds"))
nseidat          <- readRDS(here("inputs", "sablefish_data_19May2026.rds"))

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

# Recruitment deviations / age-3 abundance
ss3_age3 <- rowSums(
  SPoRC::sgl_rg_wc_sablefish_data$ss3$NAA[, 3, ]
)[ss3_year_idx]

ss3_recdev <- (ss3_age3 - mean(ss3_age3)) / sd(ss3_age3)


# SS3 is Alaska-wide, so give it a single region ------------------------
ss3_ts <- tibble(
  Region   = "WC",
  Year     = ss3_years,
  Rec      = ss3_rec,
  SSB      = ss3_ssb,
  RecDevs  = ss3_recdev
)


# NSEI model ------------------------------------------------------------
nsei_rec <- tibble(
  Region   = "CS",
  Year     = years,
  Rec      = nsei$pred_rec,
  SSB      = rowSums(nsei$spawn_biom)[-52],
  RecDevs  = (nsei$pred_rec - mean(nsei$pred_rec)) /
    sd(nsei$pred_rec)
)


# Combine ---------------------------------------------------------------
ts_df <- bind_rows(
  rtmb_rec %>%
    left_join(rtmb_ssb, by = c("Region", "Year")) %>%
    left_join(rtmb_recdev, by = c("Region", "Year")),

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
