# Loop through species across a station

ADD LATER

## Usage

``` r
loop_species_across_station(
  T_vector,
  starting_weight = 6.382417,
  prey_ED = 3698,
  oxycal_coeff = 13560
)
```

## Arguments

- T_vector:

  Temperatures (degrees C), starting at Julian day 1

- starting_weight:

  Starting weight (grams) for the species; can be a single value or a
  vector with one value per organism

- prey_ED:

  joules per gram of wet mass for the species' prey; can be a single
  value or a vector with one value per organism

- oxycal_coeff:

  joules per gram of oxygen for the species; can be a single value or a
  vector with one value per organism

- species:

  Scientific name of the species (e.g. "Micropterus salmoides")

## Value

A data.frame with growth results for each species and life stage,
including species name, life stage, and cumulative weight over the year
