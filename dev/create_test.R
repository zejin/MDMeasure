library(devtools)

# create a test workflow
if (.Platform$OS.type == "windows") {
  setwd("C:/Academia/Cornell/Research/Mutual Multivariate Independence")
} else {
  setwd("~")
}

devtools::use_testthat("MDMeasure")
