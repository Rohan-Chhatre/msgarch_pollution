###########################################################################################################################
## NY
###########################################################################################################################
ny_multi = read.csv("/Users/rohan/Documents/Rohan/Fall 25/STATINF1/data/NY_multi.csv")
ny_multi$date = as.Date(ny_multi$date)

fit_pollutant_msgarch_ny <- function(ny_multi, pollutant, k) {
  # 1. Extract raw series and dates -----------------------------------------
  x_raw <- ny_multi %>%
    arrange(date) %>%
    pull(.data[[pollutant]])
  
  x_raw <- na.omit(x_raw)
  
  dates <- ny_multi %>%
    arrange(date) %>%
    filter(!is.na(.data[[pollutant]])) %>%
    pull(date)
  
  # 2. ARIMA on levels ------------------------------------------------------
  x_ts   <- ts(x_raw, frequency = 365)
  fit_ar <- auto.arima(x_ts, seasonal = TRUE)
  x_resid <- residuals(fit_ar)
  
  # ---- 2a. ACF/PACF plots + KPSS on ARIMA residuals -----------------------
  y_diag <- x_resid
  
  # Directories for diagnostics
  acf_dir   <- "/Users/rohan/Documents/Rohan/Fall 25/STATINF1/plots/acf_residuals_NY"
  pacf_dir  <- "/Users/rohan/Documents/Rohan/Fall 25/STATINF1/plots/pacf_residuals_NY"
  kpss_dir  <- "/Users/rohan/Documents/Rohan/Fall 25/STATINF1/plots/kpss_residuals_NY"
  
  dir.create(acf_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(pacf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(kpss_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ACF plot
  acf_file <- file.path(
    acf_dir,
    paste0("ACF_residuals_", pollutant, "_NY.png")
  )
  png(acf_file, width = 1200, height = 800)
  acf(
    y_diag,
    main = paste0("ACF of ARIMA residuals - ", pollutant, " (NY)")
  )
  dev.off()
  
  # PACF plot
  pacf_file <- file.path(
    pacf_dir,
    paste0("PACF_residuals_", pollutant, "_NY.png")
  )
  png(pacf_file, width = 1200, height = 800)
  pacf(
    y_diag,
    main = paste0("PACF of ARIMA residuals - ", pollutant, " (NY)")
  )
  dev.off()
  
  # KPSS test
  kpss_result <- tseries::kpss.test(y_diag, null = "Level")
  kpss_file <- file.path(
    kpss_dir,
    paste0("KPSS_residuals_", pollutant, "_NY.txt")
  )
  sink(kpss_file)
  cat("KPSS test for ARIMA residuals -", pollutant, "(NY)\n\n")
  print(kpss_result)
  sink()
  
  # 3. Train/Test split -----------------------------------------------------
  y <- x_resid
  n <- length(y)
  train_size <- floor(0.8 * n)
  
  y_train <- y[1:train_size]
  y_test  <- y[(train_size + 1):n]
  
  dates_train <- dates[1:train_size]
  dates_test  <- dates[(train_size + 1):n]
  
  # 4. MS-GARCH -------------------------------------------------------------
  spec <- CreateSpec(
    variance.spec     = list(model = "sGARCH"),
    distribution.spec = list(distribution = "norm"),
    switch.spec       = list(do.mix = FALSE, K = k)
  )
  
  fit_ms <- FitML(spec = spec, data = y_train)
  
  # 5. Vol reconstruction ----------------------------------------------------
  vol_in_raw <- as.numeric(Volatility(fit_ms))
  
  pred <- predict(
    object = fit_ms,
    nahead = length(y_test),
    do.return.draw = FALSE
  )$vol   # numeric vector
  
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
  
  # Plot
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
      title = paste0("Predictions Based on ARIMA + MS-GARCH ", pollutant, " (NY)"),
      x = "Date",
      y = pollutant
    ) +
    theme_bw(base_size = 14) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  pred_dir <- "/Users/rohan/Documents/Rohan/Fall 25/STATINF1/plots/predictions_NY"
  dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = file.path(
      pred_dir,
      paste0("Predictions_", pollutant, "_NY.png")
    ),
    plot   = p,
    width  = 10,
    height = 6,
    dpi    = 300
  )
  
  # 7. RMSE -----------------------------------------------------------------
  rmse_insample  <- Metrics::rmse(df_plot$actual_var[1:train_size],
                                  df_plot$vol_in_var[1:train_size])
  rmse_outsample <- Metrics::rmse(df_plot$actual_var[(train_size + 1):n],
                                  df_plot$vol_out_var[(train_size + 1):n])
  
  ###########################################################################
  ## 8. NON-STATIONARY EXTREME VALUE ANALYSIS USING RtsEva
  ##    EVA on ACTUAL series (non-stationary)
  ###########################################################################
  
  timeAndSeries <- data.frame(
    date  = dates,
    value = as.numeric(x_ts)
  )
  
  timeWindow <- 10 * 365
  
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
  
  # Use last time index for GPD params
  n_time    <- length(nonStationaryEvaParams[[2]]$parameters$sigma)
  timeIndex <- n_time
  
  xi    <- nonStationaryEvaParams[[2]]$parameters$epsilon[timeIndex]
  sigma <- nonStationaryEvaParams[[2]]$parameters$sigma[timeIndex]
  thr   <- nonStationaryEvaParams[[2]]$parameters$threshold[timeIndex]
  
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
  
  rl10 <- rl10_res$ReturnLevels
  rl20 <- rl20_res$ReturnLevels
  
  eva_dir <- "/Users/rohan/Documents/Rohan/Fall 25/STATINF1/plots/eva"
  dir.create(eva_dir, recursive = TRUE, showWarnings = FALSE)
  
  eva_file <- file.path(
    eva_dir,
    paste0("EVA_RtsEva_GPD_", pollutant, "_NY.png")
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
  
  ###########################################################################
  ## RETURN -----------------------------------------------------------------
  ###########################################################################
  
  list(
    pollutant      = pollutant,
    arima_fit      = fit_ar,
    msgarch_fit    = fit_ms,
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
    
    threshold_95     = thr,
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
fit_co_ny   <- fit_pollutant_msgarch_ny(ny_multi, "CO",  k)
fit_so2_ny  <- fit_pollutant_msgarch_ny(ny_multi, "SO2", k)
fit_no2_ny  <- fit_pollutant_msgarch_ny(ny_multi, "NO2", k)
fit_pm25_ny <- fit_pollutant_msgarch_ny(ny_multi, "PM25", k)

# Plots:
fit_co_ny$plot
fit_so2_ny$plot
fit_no2_ny$plot
fit_pm25_ny$plot

# ARIMA outputs:
fit_co_ny$arima_fit
fit_so2_ny$arima_fit
fit_no2_ny$arima_fit
fit_pm25_ny$arima_fit

# MS-GARCH outputs:
fit_co_ny$msgarch_fit
fit_so2_ny$msgarch_fit
fit_no2_ny$msgarch_fit
fit_pm25_ny$msgarch_fit

# RMSE In sample
fit_co_ny$rmse_insample
fit_so2_ny$rmse_insample
fit_no2_ny$rmse_insample
fit_pm25_ny$rmse_insample

# RMSE Out sample
fit_co_ny$rmse_outsample
fit_so2_ny$rmse_outsample
fit_no2_ny$rmse_outsample
fit_pm25_ny$rmse_outsample

# Return Levels 10 Year
fit_co_ny$return_level_10
fit_so2_ny$return_level_10
fit_no2_ny$return_level_10
fit_pm25_ny$return_level_10

# Return Levels 20 Year
fit_co_ny$return_level_20
fit_so2_ny$return_level_20
fit_no2_ny$return_level_20
fit_pm25_ny$return_level_20

#
fit_co_ny$msgarch_fit
fit_so2_ny$msgarch_fit
fit_no2_ny$msgarch_fit
fit_pm25_ny$msgarch_fit