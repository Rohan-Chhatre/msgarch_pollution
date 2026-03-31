# =======================================
# CT Data Download (Daily, State-Level)
# =======================================

setwd('/Users/rohan/Documents/Rohan/Fall 25/STATINF1')

library(httr)
library(jsonlite)
library(dplyr)
library(lubridate)
library(purrr)

aqs_email <- "rohanchhatre@gmail.com"      # <-- PUT YOUR EMAIL HERE
aqs_key   <- "dunwren87"     # <-- PUT YOUR KEY HERE
get_daily_state_aqs <- function(email, key,
                                param_code,   # e.g. "42101" for CO
                                state_fips,   # e.g. "09" for CT
                                bdate, edate  # "YYYYMMDD"
) {
  base_url <- "https://aqs.epa.gov/data/api/dailyData/byState"
  
  params <- list(
    email = email,
    key   = key,
    param = param_code,
    bdate = bdate,
    edate = edate,
    state = state_fips
  )
  
  res <- httr::GET(base_url, query = params)
  httr::stop_for_status(res)
  
  txt  <- httr::content(res, as = "text", encoding = "UTF-8")
  json <- jsonlite::fromJSON(txt)
  
  if (!"Data" %in% names(json) || length(json$Data) == 0) {
    warning("No data returned for param = ", param_code,
            ", state = ", state_fips,
            ", period = ", bdate, "–", edate)
    return(NULL)
  }
  
  dplyr::as_tibble(json$Data)
}

years <- 2016:2024
ct_fips <- "09"  # Connecticut

param_codes <- c(
  CO   = "42101",
  SO2  = "42401",
  NO2  = "42602",
  O3   = "44201",
  PM25 = "88101"  # PM2.5 FRM/FEM
)

library(httr)
library(jsonlite)
library(dplyr)

get_ct_pollutant_daily <- function(poll_name, param_code,
                                   years = 2016:2024,
                                   state_fips = "09") {  # CT = "09"
  
  message("Downloading ", poll_name, " for CT (param = ", param_code, 
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
    warning("No data returned for ", poll_name, " in CT.")
    return(NULL)
  }
  
  # --- NEW: detect which column holds the measurement ---
  # Common possibilities in AQS outputs:
  #   - "sample_measurement"
  #   - "arithmetic_mean"
  #   - "obs_value" (less common, but just in case)
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


ct_co_raw   <- get_ct_pollutant_daily("CO",   param_codes["CO"])
ct_so2_raw  <- get_ct_pollutant_daily("SO2",  param_codes["SO2"])
ct_no2_raw  <- get_ct_pollutant_daily("NO2",  param_codes["NO2"])
ct_pm25_raw <- get_ct_pollutant_daily("PM25", param_codes["PM25"])

make_state_daily <- function(df, new_name) {
  df %>%
    group_by(date) %>%
    summarise(
      !!new_name := mean(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(date)
}

ct_co_daily   <- make_state_daily(ct_co_raw,   "CO")
ct_so2_daily  <- make_state_daily(ct_so2_raw,  "SO2")
ct_no2_daily  <- make_state_daily(ct_no2_raw,  "NO2")
ct_pm25_daily <- make_state_daily(ct_pm25_raw, "PM25")

ct_multi <- ct_co_daily %>%
  full_join(ct_so2_daily,  by = "date") %>%
  full_join(ct_no2_daily,  by = "date") %>%
  full_join(ct_pm25_daily, by = "date") %>%
  arrange(date)

# Path to your folder
save_dir <- "/Users/rohan/Documents/Rohan/Fall 25/STATINF1/data"

# Save CT dataset
write.csv(ct_multi,
          file = file.path(save_dir, "CT_multi.csv"),
          row.names = FALSE)
