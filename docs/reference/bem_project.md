# Project growth and performance indices across heterogeneous habitat patches

solve daily energy budgets for n habitat patches

## Usage

``` r
bem_project(
  start_M2 = 100,
  temperature,
  parms.intrinsic,
  parms.temporal,
  C_eq,
  R_eq
)
```

## Arguments

- start_M2:

  (default is 100 grams)

- temperature:

  a list of dataframes populated with a time series of mean daily water
  temperature (degrees C); each list element is a different habitat
  patch (default is 20 degrees C)

- parms.intrinsic:

  A parms.intrinsic object

- parms.temporal:

  a list of dataframes populated with a time series of intrinsic (e.g.,
  GSI) and extrinsic (e.g., prey energy density) biological parameters;
  each list element is a different habitat patch

- C_eq:

  Specify consumption equation 1, 2, or 3 from Hanson et al. 1997

- R_eq:

  Specify respiration equation 1 or 2 from Hanson et al. 1997

## Value

a list of length 2; raw bem_grow output for each patch and year-end
performance indices for each patch
