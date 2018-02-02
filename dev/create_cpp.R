library(devtools)

# set up the package with Rcpp
if (.Platform$OS.type == "windows") {
  setwd("C:/Academia/Cornell/Research/Mutual Multivariate Independence")
} else {
  setwd("~")
}

devtools::use_rcpp("MDMeasure")
