# msgarch_pollution

A time-series and extreme-value analysis project for daily air pollution data using **ARIMA**, **Markov-switching GARCH (MS-GARCH)**, and **non-stationary extreme-value methods** in R.

## Overview

This repository studies daily state-level air pollution series for four U.S. states:

- Connecticut (CT)
- Massachusetts (MA)
- New Jersey (NJ)
- New York (NY)

The workflow combines:

1. **Data acquisition and cleaning** from the EPA AQS API
2. **Exploratory time-series diagnostics**
3. **Univariate ARIMA + MS-GARCH modeling**
4. **Extreme value analysis (EVA)** using non-stationary GPD / GEV methods

The main goal is to model pollutant dynamics, capture regime-switching volatility, and estimate return levels for extreme pollution events.

---

## Repository Structure

```text
msgarch_pollution/
├── code/
│   ├── data cleaning/
│   │   ├── ct.R
│   │   ├── ma.R
│   │   ├── nj.R
│   │   └── ny.R
│   ├── data description/
│   │   └── data_description.R
│   └── univariate msgarch/
│       ├── ct_msgarch.R
│       ├── ma_msgarch.R
│       ├── nj_msgarch.R
│       └── ny_msgarch.R
├── data/
│   ├── CT_multi.csv
│   ├── MA_multi.csv
│   ├── NJ_multi.csv
│   └── NY_multi.csv
└── StatInf_1_Project (4).pdf
