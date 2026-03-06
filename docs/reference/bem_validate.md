# Project growth to validate a bioenergetics model with empirical end-of-year mass

simulate growth with resampled CP and ACT parameters to identify
parameter sets that match empirical end-of-year mass

## Usage

``` r
bem_validate(
  start_M2 = 100,
  end_M2.empirical,
  temperature = 20,
  parms.intrinsic,
  parms.temporal,
  C_eq,
  R_eq,
  resamp.n = 1000
)
```

## Arguments

- start_M2:

  (default is 100 grams)

- end_M2.empirical:

  mass of fish at end of year based on empirically-derived field or lab
  size-at-age estimate

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

- resamp.n:

  specify the number of CP-ACT parameter sets to draw from uniform
  distributions (default is 1000)

## Value

a list of length 5; raw bem_grow output for each patch and year-end
performance indices for each patch; also plotted validation output
