# =======================================
# NY Data Download (Daily, State-Level)
# =======================================

ny_fips <- "36"  # New York FIPS code

get_ny_pollutant_daily <- function(poll_name, param_code,
                                   years = 2016:2024,
                                   state_fips = ny_fips) {
  
  message("Downloading ", poll_name, " for NY (param = ", param_code, 
          ") for years: ", paste(years, collapse = ", "), " ...")
  
  out_list <- purrr::map(years, function(y) {
    bdate <- sprintf("%d0101", y)  # YYYY0101
    edate <- sprintf("%d1231", y)  # YYYY1231
    
    get_daily_state_aqs(
      email      = aqs_email,
      key        = aqs_key,
      param_code = param_code,
      state_fips = state_fips,
      bdate      = bdate,
      edate      = edate
    )
  })
  
  df <- dplyr::bind_rows(out_list)
  if (nrow(df) == 0) {
    warning("No data returned for ", poll_name, " in NY.")
    return(NULL)
  }
  
  # detect which column holds the measurement
  candidate_cols <- c("sample_measurement", "arithmetic_mean", "obs_value")
  value_col <- intersect(candidate_cols, names(df))[1]
  
  if (is.na(value_col)) {
    stop("Could not find a measurement column in the data. ",
         "Available columns are: ", paste(names(df), collapse = ", "))
  }
  
  df %>%
    dplyr::mutate(
      date  = lubridate::ymd(.data[["date_local"]]),
      value = as.numeric(.data[[value_col]])
    ) %>%
    dplyr::arrange(date)
}

# --- Download raw NY pollutant data ---
ny_co_raw   <- get_ny_pollutant_daily("CO",   param_codes["CO"])
ny_so2_raw  <- get_ny_pollutant_daily("SO2",  param_codes["SO2"])
ny_no2_raw  <- get_ny_pollutant_daily("NO2",  param_codes["NO2"])
ny_pm25_raw <- get_ny_pollutant_daily("PM25", param_codes["PM25"])

# --- Aggregate to state-level daily means ---
ny_co_daily   <- make_state_daily(ny_co_raw,   "CO")
ny_so2_daily  <- make_state_daily(ny_so2_raw,  "SO2")
ny_no2_daily  <- make_state_daily(ny_no2_raw,  "NO2")
ny_pm25_daily <- make_state_daily(ny_pm25_raw, "PM25")

# --- Combine into multivariate NY dataset ---
ny_multi <- ny_co_daily %>%
  full_join(ny_so2_daily,  by = "date") %>%
  full_join(ny_no2_daily,  by = "date") %>%
  full_join(ny_pm25_daily, by = "date") %>%
  arrange(date)

# --- Save NY dataset ---
save_dir <- "/Users/rohan/Documents/Rohan/Fall 25/STATINF1/data"

write.csv(ny_multi,
          file = file.path(save_dir, "NY_multi.csv"),
          row.names = FALSE)
rm(ny_multi)
