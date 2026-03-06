# Daily water temperature many lakes

This dataset contains modeled maximum daily water temperature for every
day of the year 2022 for all interconfluence reaches of the Little
River, TN.

## Usage

``` r
temps_lake
```

## Format

A data frame with 365 rows and 8046 columns:

- date:

  (date) dates from 1 January 2021 to 31 December 2021; note that temps
  are the mean from years 2013 to 2021

- Hylak_id:

  (numeric) mean daily water temperature in degrees C for all 8046
  unique identifier of lake, from HydroLAKES (Messager et al., 2016)

## Source

LakeTEMP; Korver et al. 2024. Remote Sensing of Environment.
https://doi.org/10.1016/j.rse.2024.114164.
