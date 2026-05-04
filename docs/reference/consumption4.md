# Consumption equation 4

Consumption equation 4; new equation for fishyenergy that simply matches
Cmax for a given temperature

## Usage

``` r
consumption4(M, T, CP = 1, parms.intrinsic, match.table)
```

## Arguments

- M:

  Mass (grams) of fish

- T:

  Temperature (degrees C) at which consumption is calculated

- CP:

  Proportion of maximum consumption (default is 1.0)

- parms.intrinsic:

  A parms.intrinsic object; note that temperature dependent parameters
  are bypassed

- match.table:

  temperature dependent parameters formatted as a table with three
  columns named WT_mean, C_match, and R_match

## Value

Specific consumption rate (grams of food consumed per gram of fish mass
per day)
