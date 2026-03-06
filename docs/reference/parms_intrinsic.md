# Create a new object with intrinsic parameters

Create a new object with intrinsic parameters

## Usage

``` r
parms_intrinsic(
  CA,
  CB,
  CTM,
  CTO,
  CQ,
  CTL,
  CK1,
  CK4,
  RA,
  RB,
  RTM,
  RTO,
  RQ,
  FA,
  UA,
  SDA,
  predED
)
```

## Arguments

- CA:

  intercept of allometric mass function

- CB:

  slope of allometric mass function; should be negative bc bigger fish
  consume less per gram of body mass than smaller fish

- CTM:

  critical thermal maximum (degrees C)

- CTO:

  laboratory temperature preferendum (degrees C)

- CQ:

  approximates a Q10; the rate at which the function increases over
  relatively low water temperatures

- CTL:

  temperature-dependent parameter for consumption

- CK1:

  temperature-dependent parameter for consumption

- CK4:

  temperature-dependent parameter for consumption

- RA:

  intercept of allometric activity function

- RB:

  slope of allometric activity function

- RTM:

  critical thermal maximum (degrees C)

- RTO:

  laboratory temperature preferendum (degrees C)

- RQ:

  approximates a Q10; the rate at which the function increases over
  relatively low water temperatures

- FA:

  fish body mass (grams)

- UA:

  activity (mg/kg-day)

- SDA:

  specific dynamic action (proportion)

- predED:

  Energy density of the predator (i.e., focal) species in joules per
  gram of wet mass

## Value

An intrinsic parameters object
