library(devtools)

# check the downstream dependencies
if (.Platform$OS.type == "windows") {
  setwd("C:/Academia/Cornell/Research/Mutual Multivariate Independence/Simulation")
} else {
  setwd("~")
}

devtools::revdep_check("MDMeasure")