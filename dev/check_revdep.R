library(devtools)

# check the downstream dependencies
if (.Platform$OS.type == "windows") {
  setwd("C:/Academia/Cornell/Research/Mutual Multivariate Independence")
} else {
  setwd("~")
}

devtools::revdep_check("MDMeasure")