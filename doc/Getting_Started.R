## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(fishyenergy)

## -----------------------------------------------------------------------------
bem_raw <- read.csv(system.file("extdata","BEM_parameters_v01.csv", package="fishyenergy"))

