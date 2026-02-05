# Getting Started

The model here is adopted from <https://github.com/troiamj/GB_BEM_v01>,
which is based on [Troia and Perkin
(2022)](https://doi.org/10.1093/conphys/coac035).

``` r

library(fishyenergy)
```

First, let’s load in our parameters for our model. For the vignette,
we’re loading it from a saved one in our package, but you can get one
from your computer instead

``` r

bem_raw <- read.csv(system.file("extdata","BEM_parameters_v01.csv", package="fishyenergy"))
chosen_row <- 1
bem_vector <- bem_raw[chosen_row,]
BEM <- bem_new(bem_vector$CP, bem_vector$CA, bem_vector$CB, bem_vector$CTM, bem_vector$CTO, bem_vector$CQ, bem_vector$ACT, bem_vector$RA, bem_vector$RB, bem_vector$RTM, bem_vector$RTO, bem_vector$RQ, bem_vector$FA, bem_vector$UA, bem_vector$SDA, bem_vector$ED)
```

And then load in the temperatures for a year at a site:

``` r

temp_raw <- read.csv(system.file("extdata","modeled_WT_HUC12100201.csv", package="fishyenergy"))
chosen_col <- 1
temp_vec <- temp_raw[,chosen_col]
```

Let’s plot the temperature data to make sure they look like what we’d
expect

``` r

plot(x=as.Date(0:364, origin="2020-01-01"), y=temp_vec, xlab="Day", ylab="Temperature", type='l', bty="n")
```

![Line plot showing temperature over a
year](Getting_Started_files/figure-html/unnamed-chunk-4-1.png)

Yay, looks good. Now let’s compute the weight gain per day of our fishy:

``` r

station_result <- single_station_compute(T_vector=temp_vec, BEM=BEM)
```

And look at how its weight changes over the year at that station:

``` r

plot(x=as.Date(0:364, origin="2020-01-01"), y=station_result$W2_cum, xlab="Day", ylab="Mass (g)", type='l', bty="n")
```

![Line plot showing mass over a
year](Getting_Started_files/figure-html/unnamed-chunk-6-1.png)
