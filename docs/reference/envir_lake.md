# Daily water temperature many lakes

This dataset contains modeled maximum daily water temperature for every
day of the year 2022 for all interconfluence reaches of the Little
River, TN.

## Usage

``` r
envir_lake
```

## Format

A data frame with 365 rows and 8046 columns:

- Hylak_id:

  (numeric) mean daily water temperature in degrees C for all 8046
  unique identifier of lake, from HydroLAKES (Messager et al., 2016)

- center_long:

  (numeric) longitude of lake center in decimal degrees

- center_lat:

  (numeric) latitude of lake center in decimal degrees

- stat_method:

  (numeric) see Korver et al. 2024 for description

- n_obs:

  (numeric) see Korver et al. 2024 for description

- intermittency:

  (numeric) see Korver et al. 2024 for description

- Lake_type:

  (numeric) 1 = natural lake, 2 = reservoir, 3 = natural lake with
  regulation structure

- Elevation:

  (numeric) lake surface elevation in meters above sea level

- Lake_area:

  (numeric) lake surface area in square kilometers

- Depth_avg:

  (numeric) mean lake depth in meters

## Source

LakeTEMP; Korver et al. 2024. Remote Sensing of Environment.
https://doi.org/10.1016/j.rse.2024.114164.
