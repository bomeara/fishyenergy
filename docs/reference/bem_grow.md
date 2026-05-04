# Solve daily energy budgets to simulate growth over a calendar year

the Wisconsin bioenergetics model from Hanson et al. 1997

## Usage

``` r
bem_grow(
  start_M2 = 100,
  temperature,
  parms.intrinsic,
  parms.temporal,
  C_eq,
  R_eq,
  match.table
)
```

## Arguments

- start_M2:

  (default is 100 grams)

- temperature:

  a dataframe populated with a time series of mean daily water
  temperature (degrees C x 10) of a habitat patch

- parms.intrinsic:

  A parms.intrinsic object; note that temperature dependent parameters
  are bypassed if C_eq = 4 and/or R_eq = 4

- parms.temporal:

  a dataframe populated with a time series of intrinsic (e.g., GSI) and
  extrinsic (e.g., prey energy density) biological parameters

- C_eq:

  Specify consumption equation 1, 2, or 3 from Hanson et al. 1997 or
  equation 4

- R_eq:

  Specify respiration equation 1 or 2 from Hanson et al. 1997

- match.table:

  temperature dependent parameters formatted as a table with three
  columns named WT_mean, C_match, and R_match; only applicable if C_eq =
  4 and/or R_eq = 4

## Value

dataframe populated with output parameters (columns) for time series of
projected days (rows)
