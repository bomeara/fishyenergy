# Compute growth for a year at a station

ADD LATER

## Usage

``` r
single_station_compute(
  T_vector,
  BEM,
  starting_weight = 6.382417,
  prey_ED = 3698,
  oxycal_coeff = 13560
)
```

## Arguments

- T_vector:

  Temperatures (degrees C), starting at Julian day 1

- BEM:

  BEM object

- starting_weight:

  Starting weight (grams) from Groeschel-Taylor et al

- prey_ED:

  joules per gram of wet mass

- oxycal_coeff:

  joules per gram of oxygen

## Value

data.frame with growth over year
