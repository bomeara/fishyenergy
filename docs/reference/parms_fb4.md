# Intrinsic bioenergetics parameters from FB4

This dataset contains intrinsic bioenergetic parameters for five
species.

## Usage

``` r
parms_fb4
```

## Format

A data frame with 5 rows and 27 columns:

- genus_species:

  (character) description

- LifeStage:

  (character) description

- Source:

  (character) description

- lwA:

  (numeric) y-intercept for length-weight relationship on log-log scale

- lwB:

  (numeric) slope for length-weight relationship on log-log scale

- CEQ:

  (integer) consumption equation 1, 2, 3, or 4

- CA:

  (numeric) intercept estimate for mass-dependent consumption

- CB:

  (numeric) slope estimate for mass-dependent consumption

- CQ:

  (numeric) Q10 for temperature-dependent consumption

- CTO:

  (numeric) thermal optimum for consumption

- CTM:

  (numeric) thermal maximum for consumption

- CTL:

  (numeric) temperature-dependent parameter for consumption

- CK1:

  (numeric) temperature-dependent parameter for consumption

- CK4:

  (numeric) temperature-dependent parameter for consumption

- REQ:

  (integer) respiration equation 1, 2, or 4

- RA:

  (numeric) intercept estimate for mass-dependent respiration

- RB:

  (numeric) slope estimate for mass-dependent respiration

- RQ:

  (numeric) Q10 for temperature-dependent respiration

- RTO:

  (numeric) thermal optimum for respiration

- RTM:

  (numeric) thermal maximum for respiration

- RTL:

  (numeric) temperature-dependent parameter for respiration

- RK1:

  (numeric) temperature-dependent parameter for respiration

- RK4:

  (numeric) temperature-dependent parameter for respiration

- SDA:

  (numeric) specific dynamic action

- FA:

  (numeric) egestion

- UA:

  (numeric) excretion

- pred_ED:

  (numeric) energy density of predator (i.e., focal species)

## Source

Fish Bioenergetics 4 (Deslauriers et al. 2017).
