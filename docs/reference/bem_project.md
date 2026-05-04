# Project performance indices across heterogeneous habitat patches

solve daily energy budgets for n habitat patches

## Usage

``` r
bem_project(
  start_M2 = 100,
  temperature,
  parms.intrinsic,
  parms.temporal,
  C_eq,
  R_eq,
  is.raster = FALSE
)
```

## Arguments

- start_M2:

  (default is 100 grams)

- temperature:

  a dataframe with 365 rows; the first column is the date; subsequent
  columns are WT_mean (degrees C x 10) for n habitat patches OR a raster
  layer with x, y, and z values representing latitudes, longitudes, and
  365 days of temperature

- parms.intrinsic:

  A parms.intrinsic object

- parms.temporal:

  a list of dataframes populated with a time series of intrinsic (e.g.,
  GSI) and extrinsic (e.g., prey energy density) biological parameters;
  each list element is a different habitat patch

- C_eq:

  Specify consumption equation 1, 2, or 3 from Hanson et al. 1997 or
  equation 4

- R_eq:

  Specify respiration equation 1 or 2 from Hanson et al. 1997

- is.raster:

  FALSE if temperature data are in tabular format; TRUE if temperature
  data are in raster format. Default is FALSE

## Value

a list of length 2; raw bem_grow output for each patch and year-end
performance indices for each patch
