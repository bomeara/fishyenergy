# Respiration equation 4

Respiration equation 4; new equation for fishyenergy that simply matches
Rrest for a given temperature

## Usage

``` r
respiration4(M, T, ACT = 1, parms.intrinsic, match.table)
```

## Arguments

- M:

  Mass (grams) of fish

- T:

  Temperature (degrees C) at which consumption is calculated

- ACT:

  Activity multiplier (default is 1.0)

- parms.intrinsic:

  A parms.intrinsic object; note that temperature dependent parameters
  are bypassed

- match.table:

  temperature dependent parameters formatted as a table with three
  columns named WT_mean, C_match, and R_match

## Value

Specific respiration rate (grams of oxygen per gram of fish mass per
day)
