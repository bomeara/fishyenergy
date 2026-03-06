# Create a list object of temperature time series for habitat patches

reformats a data frame to a list for use in bem_project function

## Usage

``` r
compile_temps(temperature, is.raster = FALSE)
```

## Arguments

- temperature:

  a dataframe with 365 rows; the first column is the date; subsequent
  columns are WT_mean for n habitat patches

## Value

a list of dataframes populated with a time series of mean daily water
temperature (degrees C); each list element is a different habitat patch
(default is 20 degrees C)
