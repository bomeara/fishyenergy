# Create a BEM object from FB4 data

Create a BEM object from FB4 data

## Usage

``` r
bem_from_fb4(
  species,
  lifestage = NULL,
  CP = 0.5,
  FA = 0.1018428,
  UA = 0.08743042,
  SDA = 0.150819,
  ED = 3983.646
)
```

## Arguments

- species:

  Scientific name of the species (e.g. "Micropterus salmoides")

- lifestage:

  Life stage of the species (e.g. "Adult"). If NA, uses the first life
  stage found for the species.

- CP:

  proportion of maximum consumption (default is 0.5)

- FA:

  egestion proportion (what proportion of food gets passed through)
  (default is 0.1018428)

- UA:

  excretion proportion (what proportion is assimilated but not allocated
  towards metabolism or gonad growth) (default is 0.08743042)

- SDA:

  specific dynamic action, proportion of energy going to digestion
  (default is 0.150819)

- ED:

  energy density (J/kg) (default is 3983.646)

## Value

BEM object with parameters for the specified species and life stage
