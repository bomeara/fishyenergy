# Environmental covariates for many stream reaches

This dataset contains geographic and hydrographic covariates for all 513
interconfluence stream reaches in the Little River, TN.

## Usage

``` r
envir_casc
```

## Format

A vector shapefile with an attribute table containing 513 rows and 5
columns:

- COMID:

  (character) Common Identifier for NHD reaches

- GNIS_NAME:

  (character) Stream name from NHD after Geographic Names Information
  System

- FCode:

  (integer) Flowline type from NHD

- TotDASqKM:

  (numeric) Catchment area draining to the interconfluence stream reach
  in square kilometers

- ElevCat:

  (numeric) Median elevation of the interconfluence stream reach in
  meters above sea level

## Source

Southeastern Climate Adaptation Science Center (CASC), National
Hydrography Dataset (NHD), StreamCat
