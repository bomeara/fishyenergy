# Plot temperature dependent curves

Plot the thermal reaction norms for consumption, respiration, and other
intrinsic rates

## Usage

``` r
bem_curve(M, T, CP = 1, ACT = 1, prey_ED = 4000, parms.intrinsic, C_eq, R_eq)
```

## Arguments

- M:

  Mass (grams) of fish

- T:

  Vector of temperatures representing the temperature range to be
  plotted

- CP:

  Proportion of maximum consumption (default is 1.0)

- ACT:

  Activity multiplier (default is 1.0)

- prey_ED:

  Prey energy density in joules per gram of wet mass (default is 4000)

- parms.intrinsic:

  A parms.intrinsic object

- C_eq:

  Specify consumption equation 1, 2, or 3 from Hanson et al. 1997

- R_eq:

  Specify respiration equation 1 or 2 from Hanson et al. 1997

## Value

A plot of temperature dependent rates
