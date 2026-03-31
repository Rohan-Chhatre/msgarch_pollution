###########################################################################################################################
##CT
###########################################################################################################################

library(tidyverse)
library(forecast)
library(MSGARCH)
library(ggplot2)
library(tseries)
library(RtsEva)

ct_multi = read.csv("/Users/rohan/Documents/Rohan/Fall 25/STATINF1/data/CT_multi.csv")
ct_multi$date = as.Date(ct_multi$date)

fit_pollutant_msgarch <- function(ct_multi, pollutant, k) {
  # 1. Extract raw series and dates -----------------------------------------
  x_raw <- ct_multi %>%
    arrange(date) %>%
    pull(.data[[pollutant]])
  
  x_raw <- na.omit(x_raw)
  
  dates <- ct_multi %>%
    arrange(date) %>%
    filter(!is.na(.data[[pollutant]])) %>%
    pull(date)
  
  # 2. ARIMA on levels ------------------------------------------------------
  x_ts   <- ts(x_raw, frequency = 365)
  fit_ar <- auto.arima(x_ts, seasonal = TRUE)
  x_resid <- residuals(fit_ar)
  
  # ---- 2a. ACF/PACF plots + KPSS on ARIMA residuals -----------------------
  # Use full residual series y for these diagnostics
  y_diag <- x_resid
  
  # Directories for diagnostics
  acf_dir   <- "/Users/rohan/Documents/Rohan/Fall 25/STATINF1/plots/acf_residuals"
  pacf_dir  <- "/Users/rohan/Documents/Rohan/Fall 25/STATINF1/plots/pacf_residuals"
  kpss_dir  <- "/Users/rohan/Documents/Rohan/Fall 25/STATINF1/plots/kpss_residuals"
  
  dir.create(acf_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(pacf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(kpss_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ACF plot (base R)
  acf_file <- file.path(
    acf_dir,
    paste0("ACF_residuals_", pollutant, "_CT.png")
  )
  png(acf_file, width = 1200, height = 800)
  acf(
    y_diag,
    main = paste0("ACF of ARIMA residuals - ", pollutant, " (CT)")
  )
  dev.off()
  
  # PACF plot (base R)
  pacf_file <- file.path(
    pacf_dir,
    paste0("PACF_residuals_", pollutant, "_CT.png")
  )
  png(pacf_file, width = 1200, height = 800)
  pacf(
    y_diag,
    main = paste0("PACF of ARIMA residuals - ", pollutant, " (CT)")
  )
  dev.off()
  
  # KPSS test on residuals and save output
  kpss_result <- tseries::kpss.test(y_diag, null = "Level")
  kpss_file <- file.path(
    kpss_dir,
    paste0("KPSS_residuals_", pollutant, "_CT.txt")
  )
  sink(kpss_file)
  cat("KPSS test for ARIMA residuals -", pollutant, "(CT)\n\n")
  print(kpss_result)
  sink()
  
  # 3. Split residuals into train / test -----------------------------------
  y <- x_resid
  n <- length(y)
  train_size <- floor(0.8 * n)
  
  y_train <- y[1:train_size]
  y_test  <- y[(train_size + 1):n]
  
  dates_train <- dates[1:train_size]
  dates_test  <- dates[(train_size + 1):n]
  
  # 4. MS-GARCH on train residuals -----------------------------------------
  spec <- CreateSpec(
    variance.spec     = list(model = "sGARCH"),
    distribution.spec = list(distribution = "norm"),
    switch.spec       = list(do.mix = FALSE, K = k)
  )
  
  fit_ms <- FitML(spec = spec, data = y_train)
  
  # 5. In-sample volatility + forecast volatility --------------------------
  # In-sample sigma_t
  vol_in_raw <- as.numeric(Volatility(fit_ms))
  
  # Forecast sigma_t for test (no draws -> vector)
  pred <- predict(
    object = fit_ms,
    nahead = length(y_test),
    do.return.draw = FALSE
  )$vol   # this is already a numeric vector
  
  # Your reconstruction: sign(y) * sigma + ARIMA fitted mean
  vol_in <- sign(y_train) * vol_in_raw +
    fitted(fit_ar)[1:train_size]
  
  vol_forecast <- sign(y_test) * pred +
    fitted(fit_ar)[(train_size + 1):n]
  
  # 6. Build plot data ------------------------------------------------------
  df_plot <- tibble(
    date        = dates,
    actual_var  = x_ts,
    vol_in_var  = c(vol_in,           rep(NA, n - train_size)),
    vol_out_var = c(rep(NA, train_size), vol_forecast)
  )
  df_plot$date = as.Date(df_plot$date)
  
  # --- 6a. GGplot with legend & save to disk -------------------------------
  p <- ggplot(df_plot, aes(x = date)) +
    geom_line(aes(y = actual_var,  colour = "Actual"), alpha = 0.6) +
    geom_line(aes(y = vol_in_var,  colour = "In-sample prediction"),
              size = 0.8, na.rm = TRUE) +
    geom_line(aes(y = vol_out_var, colour = "Out-of-sample prediction"),
              size = 0.9, na.rm = TRUE) +
    scale_colour_manual(
      name   = "Series",
      values = c(
        "Actual"                   = "grey60",
        "In-sample prediction"     = "blue",
        "Out-of-sample prediction" = "red"
      )
    ) +
    labs(
      title = paste0("Predictions Based on ARIMA + MS-GARCH ", pollutant, " (CT)"),
      x = "Date",
      y = pollutant
    ) +
    theme_bw(base_size = 14) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  # Save prediction plot
  pred_dir <- "/Users/rohan/Documents/Rohan/Fall 25/STATINF1/plots/predictions"
  dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = file.path(
      pred_dir,
      paste0("Predictions_", pollutant, "_CT.png")
    ),
    plot   = p,
    width  = 10,
    height = 6,
    dpi    = 300
  )
  
  # 7. RMSE
  rmse_insample = Metrics::rmse(df_plot$actual_var[1:train_size],
                                df_plot$vol_in_var[1:train_size])
  rmse_outsample = Metrics::rmse(df_plot$actual_var[(1+train_size):n],
                                 df_plot$vol_out_var[(1+train_size):n])
  
  ###############################################################
  ## 7. EXTREME VALUE ANALYSIS (non-stationary GPD via RtsEva)
  ##    NOTE: EVA is now done on the ACTUAL series using TsEvaNs
  ###############################################################
  
  # Build time series for RtsEva (actual pollutant series)
  timeAndSeries <- data.frame(
    date  = dates,
    value = as.numeric(x_ts)
  )
  
  # Choose a time window (e.g., 10-year window for trend / scale)
  timeWindow <- 10 * 365
  
  # Non-stationary EVA using transformed-stationary approach
  # (you must have library(RtsEva) loaded)
  tseva <- TsEvaNs(
    timeAndSeries           = timeAndSeries,
    timeWindow              = timeWindow,
    transfType              = "trendPeaks",
    minPeakDistanceInDays   = 10,
    seasonalityVar          = NA,
    minEventsPerYear        = -1,
    gevMaxima               = "annual",
    ciPercentile            = 90,
    gevType                 = "GEV",
    evdType                 = c("GEV", "GPD"),
    tail                    = "high",
    epy                     = -1,
    lowdt                   = 7,
    trans                   = "ori"
  )
  
  nonStationaryEvaParams  <- tseva[[1]]
  stationaryTransformData <- tseva[[2]]
  
  # GPD parameters are usually in the 2nd element of nonStationaryEvaParams
  # Use the last time index (most recent period) for return levels
  n_time    <- length(nonStationaryEvaParams[[2]]$parameters$sigma)
  timeIndex <- n_time
  
  # Extract GPD parameters at this time index
  xi    <- nonStationaryEvaParams[[2]]$parameters$epsilon[timeIndex]
  sigma <- nonStationaryEvaParams[[2]]$parameters$sigma[timeIndex]
  thr   <- nonStationaryEvaParams[[2]]$parameters$threshold[timeIndex]
  
  # Compute 10- and 20-year return levels (GEV+GPD) at this time index
  rl10_res <- tsEvaComputeRLsGEVGPD(
    nonStationaryEvaParams = nonStationaryEvaParams,
    RPgoal                 = 10,
    timeIndex              = timeIndex,
    trans                  = "ori"
  )
  
  rl20_res <- tsEvaComputeRLsGEVGPD(
    nonStationaryEvaParams = nonStationaryEvaParams,
    RPgoal                 = 20,
    timeIndex              = timeIndex,
    trans                  = "ori"
  )
  
  # Take the GPD-based return levels
  rl10 <- rl10_res$ReturnLevels
  rl20 <- rl20_res$ReturnLevels
  
  # Optional: save a non-stationary RL plot for GPD
  eva_dir <- "/Users/rohan/Documents/Rohan/Fall 25/STATINF1/plots/eva"
  dir.create(eva_dir, recursive = TRUE, showWarnings = FALSE)
  
  eva_file <- file.path(
    eva_dir,
    paste0("EVA_RtsEva_GPD_", pollutant, "_CT.png")
  )
  
  png(eva_file, width = 1200, height = 900)
  plt <- tsEvaPlotReturnLevelsGPDFromAnalysisObj(
    nonStationaryEvaParams  = nonStationaryEvaParams,
    stationaryTransformData = stationaryTransformData,
    timeIndex               = timeIndex,
    trans                   = "ori"
  )
  print(plt)
  dev.off()
  
  # 8. Return everything you might want ------------------------------------
  list(
    pollutant      = pollutant,
    arima_fit      = fit_ar,   # ARIMA output
    msgarch_fit    = fit_ms,   # MS-GARCH output
    df_plot        = df_plot,
    plot           = p,
    train_size     = train_size,
    dates_train    = dates_train,
    dates_test     = dates_test,
    y_train        = y_train,
    y_test         = y_test,
    vol_in         = vol_in,
    vol_forecast   = vol_forecast,
    rmse_insample  = rmse_insample,
    rmse_outsample = rmse_outsample,
    kpss_result    = kpss_result,  # KPSS on ARIMA residuals
    
    # EVA outputs (now from RtsEva, non-stationary GPD on actuals)
    threshold_95     = thr,     # effective GPD threshold at timeIndex
    gpd_fit          = list(
      nonStationaryEvaParams  = nonStationaryEvaParams,
      stationaryTransformData = stationaryTransformData,
      timeIndex               = timeIndex
    ),
    gpd_shape        = xi,
    gpd_scale        = sigma,
    return_level_10  = rl10,
    return_level_20  = rl20
  )
}

# assume k already defined, e.g. 
k <- 3
fit_co  <- fit_pollutant_msgarch(ct_multi, "CO",  k)
fit_so2 <- fit_pollutant_msgarch(ct_multi, "SO2", k)
fit_no2 <- fit_pollutant_msgarch(ct_multi, "NO2", k)
fit_pm25 <- fit_pollutant_msgarch(ct_multi, "PM25", k)

# Plots:
fit_co$plot
fit_so2$plot
fit_no2$plot
fit_pm25$plot

# ARIMA outputs:
fit_co$arima_fit
fit_so2$arima_fit
fit_no2$arima_fit
fit_pm25$arima_fit

# MS-GARCH outputs:
fit_co$msgarch_fit
fit_so2$msgarch_fit
fit_no2$msgarch_fit
fit_pm25$msgarch_fit

# RMSE In sample
fit_co$rmse_insample
fit_so2$rmse_insample
fit_no2$rmse_insample
fit_pm25$rmse_insample

# RMSE Out sample
fit_co$rmse_outsample
fit_so2$rmse_outsample
fit_no2$rmse_outsample
fit_pm25$rmse_outsample


# Return Levels 10 Year
fit_co$return_level_10
fit_so2$return_level_10
fit_no2$return_level_10
fit_pm25$return_level_10

# Return Levels 20 Year
fit_co$return_level_20
fit_so2$return_level_20
fit_no2$return_level_20
fit_pm25$return_level_20

# 
fit_co$msgarch_fit
fit_so2$msgarch_fit
fit_no2$msgarch_fit
fit_pm25$msgarch_fit
