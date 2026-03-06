# Solve daily energy budgets to simulate growth over a calendar year

the Wisconsin bioenergetics model from Hanson et al. 1997

## Usage

``` r
bem_grow(
  start_M2 = 100,
  temperature = 20,
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

  a dataframe populated with a time series of mean daily water
  temperature (degrees C) of a habitat patch (default is 20 degrees C)

- parms.intrinsic:

  A parms.intrinsic object

- parms.temporal:

  a dataframe populated with a time series of intrinsic (e.g., GSI) and
  extrinsic (e.g., prey energy density) biological parameters

- C_eq:

  Specify consumption equation 1, 2, or 3 from Hanson et al. 1997

- R_eq:

  Specify respiration equation 1 or 2 from Hanson et al. 1997

## Value

dataframe populated with output parameters (columns) for time series of
projected days (rows)
