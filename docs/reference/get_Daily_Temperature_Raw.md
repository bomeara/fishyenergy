# Get daily temperature data from monitoring stations

This uses the dataRetrieval package to fetch daily temperature data

## Usage

``` r
daily_temperature_raw(site_code, start_date, end_date)
```

## Arguments

- site_code:

  WQP site code (e.g. "USGS-12345678")

- start_date:

  Start date for the data retrieval (e.g. "2020-01-01")

- end_date:

  End date for the data retrieval (e.g. "2020-12-31")

## Value

A data frame with daily temperature data
