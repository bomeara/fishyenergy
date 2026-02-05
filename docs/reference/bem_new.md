# Create a new object with BEM data

Create a new object with BEM data

## Usage

``` r
bem_new(
  CP,
  CA,
  CB,
  CTM,
  CTO,
  CQ,
  CTL,
  CK1,
  CK4,
  ACT,
  RA,
  RB,
  RTM,
  RTO,
  RQ,
  FA,
  UA = 0.1,
  SDA = 0.15,
  ED = 3000
)
```

## Arguments

- CP:

  proportion of maximum consumption

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

  temperature at which consumption is XXXX (degrees C)

- CK1:

  temperature at which consumption is XXXX (degrees C)

- CK4:

  temperature at which consumption is XXXX (degrees C)

- ACT:

  activity coefficient

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

- ED:

  energy density (J/kg)

## Value

BEM object
