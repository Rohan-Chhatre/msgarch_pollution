library(tidyverse)
library(lubridate)
library(tseries)
library(strucchange)

diag_pollutant_ts <- function(data,
                              pollutant,
                              state,                # ← ADD ONLY THIS ARGUMENT
                              date_col   = "date",
                              freq       = 365,
                              max_breaks = 5,
                              min_seg    = 0.10,
                              do_plots   = TRUE) {
  
  # 1. Extract and clean series --------------------------------------------
  df <- data %>%
    arrange(.data[[date_col]]) %>%
    select(date = all_of(date_col), value = all_of(pollutant)) %>%
    filter(!is.na(value))
  
  if (nrow(df) == 0) {
    stop("No non-missing values for pollutant: ", pollutant)
  }
  
  n    <- nrow(df)
  y    <- df$value
  y_ts <- ts(y, frequency = freq)
  y2   <- y^2  # for variance breaks
  
  # 2. KPSS test ------------------------------------------------------------
  kpss_res <- tseries::kpss.test(y_ts, null = "Level")
  
  # 3. Bai–Perron test for variance breaks ---------------------------------
  h_val <- floor(min_seg * n)
  if (h_val < 5L) h_val <- 5L
  
  bpvar_full <- strucchange::breakpoints(y2 ~ 1, h = h_val, breaks = max_breaks)
  
  bic_vals_var <- BIC(bpvar_full)
  best_k_var   <- which.min(bic_vals_var) - 1L
  
  bpvar_final   <- breakpoints(bpvar_full, breaks = best_k_var)
  break_inds    <- bpvar_final$breakpoints
  break_dates   <- if (!all(is.na(break_inds))) df$date[break_inds] else NULL
  
  # 4. Plots & Saving -------------------------------------------------------
  if (do_plots) {
    
    ### ---- 4A. Time Series Plot ---- ###
    png(
      filename = sprintf("/Users/rohan/Documents/Rohan/Fall 25/STATINF1/plots/actual/ts_%s_%s.png",
                         state, pollutant),
      width = 600, height = 400
    )
    title = paste0("Time Series for ", pollutant, " (", state, ")")
    plot(as.Date(df$date), y_ts, col='black', ylab=pollutant, xlab='Date',
         main=title, type='l', lwd=1)
    grid(col='lightgray')
    dev.off()
    
    
    ### ---- 4B. ACF Plot ---- ###
    png(
      filename = sprintf("/Users/rohan/Documents/Rohan/Fall 25/STATINF1/plots/acf/acf_%s_%s.png",
                         state, pollutant),
      width = 600, height = 400
    )
    acf(y_ts,
        main = paste("ACF of", pollutant, "(", state, ")"),
        na.action = na.pass)
    dev.off()
    
    
    ### ---- 4C. PACF Plot ---- ###
    png(
      filename = sprintf("/Users/rohan/Documents/Rohan/Fall 25/STATINF1/plots/pacf/pacf_%s_%s.png",
                         state, pollutant),
      width = 600, height = 400
    )
    pacf(y_ts,
         main = paste("PACF of", pollutant, "(", state, ")"),
         na.action = na.pass)
    dev.off()
    
    
    ### ---- 4D. Bai–Perron Plot ---- ###
    png(
      filename = sprintf("/Users/rohan/Documents/Rohan/Fall 25/STATINF1/plots/bp/bp_%s_%s.png",
                         state, pollutant),
      width = 600, height = 400
    )
    plot(bpvar_full,
         main = paste("Bai–Perron Variance Breaks for", pollutant, "(", state, ")"))
    dev.off()
    
    
    ### ---- 4E. Variance Series with Regimes ---- ###
    png(
      filename = sprintf("/Users/rohan/Documents/Rohan/Fall 25/STATINF1/plots/bp/bp_fit_%s_%s.png",
                         state, pollutant),
      width = 600, height = 400
    )
    plot(ts(y2, frequency = freq),
         main = paste(pollutant, "^2 with Fitted Variance Regimes (", state, ")", sep=""),
         ylab = paste0(pollutant, "^2"),
         xlab = "Time")
    lines(fitted(bpvar_final), col = "red", lwd = 2)
    if (!is.null(break_inds)) {
      abline(v = break_inds, col = "blue", lty = 2)
    }
    dev.off()
    
  }
  
  # 5. Return diagnostics ---------------------------------------------------
  list(
    pollutant        = pollutant,
    n_obs            = n,
    kpss_result      = kpss_res,
    breakpoints_var_full  = bpvar_full,
    breakpoints_var_final = bpvar_final,
    break_indices_var     = break_inds,
    break_dates_var       = break_dates
  )
}

# For CT
ct_co_diag   <- diag_pollutant_ts(ct_multi, "CO", state='CT')
ct_no2_diag  <- diag_pollutant_ts(ct_multi, "NO2", state='CT')
ct_pm25_diag <- diag_pollutant_ts(ct_multi, "PM25", state='CT')
ct_so2_diag <- diag_pollutant_ts(ct_multi, "SO2", state='CT')

ct_co_diag$kpss_result
ct_no2_diag$kpss_result
ct_pm25_diag$kpss_result
ct_so2_diag$kpss_result

# For NY
ny_co_diag   <- diag_pollutant_ts(ny_multi, "CO", state='NY')
ny_no2_diag  <- diag_pollutant_ts(ny_multi, "NO2", state='NY')
ny_pm25_diag <- diag_pollutant_ts(ny_multi, "PM25", state='NY')
ny_so2_diag <- diag_pollutant_ts(ny_multi, "SO2", state='NY')

ny_co_diag$kpss_result
ny_no2_diag$kpss_result
ny_pm25_diag$kpss_result
ny_so2_diag$kpss_result

# For NJ
nj_co_diag   <- diag_pollutant_ts(nj_multi, "CO", state='NJ')
nj_no2_diag  <- diag_pollutant_ts(nj_multi, "NO2", state='NJ')
nj_pm25_diag <- diag_pollutant_ts(nj_multi, "PM25", state='NJ')
nj_so2_diag <- diag_pollutant_ts(nj_multi, "SO2", state='NJ')

nj_co_diag$kpss_result
nj_no2_diag$kpss_result
nj_pm25_diag$kpss_result
nj_so2_diag$kpss_result

# For MA
ma_co_diag   <- diag_pollutant_ts(ma_multi, "CO", state='MA')
ma_no2_diag  <- diag_pollutant_ts(ma_multi, "NO2", state='MA')
ma_pm25_diag <- diag_pollutant_ts(ma_multi, "PM25", state='MA')
ma_so2_diag <- diag_pollutant_ts(ma_multi, "SO2", state='MA')

ma_co_diag$kpss_result
ma_no2_diag$kpss_result
ma_pm25_diag$kpss_result
ma_so2_diag$kpss_result

