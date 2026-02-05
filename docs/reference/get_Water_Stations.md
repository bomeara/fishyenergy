# Get water quality monitoring stations

This can retrieve water quality monitoring stations based on state or
bounding box coordinates. It uses the dataRetrieval package to fetch the
data.

## Usage

``` r
water_stations_find(
  state = NULL,
  min_lon = NULL,
  max_lon = NULL,
  min_lat = NULL,
  max_lat = NULL
)
```

## Arguments

- state:

  State code (e.g. "CA" for California) to filter stations by state

- min_lon:

  Minimum longitude for bounding box (optional)

- max_lon:

  Maximum longitude for bounding box (optional)

- min_lat:

  Minimum latitude for bounding box (optional)

- max_lat:

  Maximum latitude for bounding box (optional) \#' @return A data frame
  with water quality monitoring stations, including site number, site
  name, latitude, longitude, and other relevant information
