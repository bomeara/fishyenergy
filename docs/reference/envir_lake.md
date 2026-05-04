# Environmental covariates for many lakes

This dataset contains geographic and bathymetric covariates for all
8,046 lakes extracted from the lakeTEMP dataset.

## Usage

``` r
envir_lake
```

## Format

A data frame with 8,046 rows and 10 columns:

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
