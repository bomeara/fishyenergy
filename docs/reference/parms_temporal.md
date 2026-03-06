# Temporally varying parameters

This dataset contains bioenergetics parameters that potentially vary
each day of a 365-day calendar year.

## Usage

``` r
parms_temporal
```

## Format

A data frame with 365 rows and 5 columns:

- date:

  (date) dates from 1 to 365

- CP:

  (numeric) Proportion of consumption; set to 1.0 assumes prey is
  unlimited

- ACT:

  (numeric) Activity multiplier; set to 1.0 assumes predator (i.e.,
  focal fish) is at rest

- gsi_f:

  (numeric) Gonadosomatic index for females; set to 0 assumes
  reproductive immaturity

- gsi_m:

  (numeric) Gonadosomatic index for males; set to 0 assumes reproductive
  immaturity

## Source

n/a.
