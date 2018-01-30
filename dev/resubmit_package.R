library(devtools)

# resubmit the source package to cran 
if (.Platform$OS.type == "windows") {
  setwd("C:/Academia/Cornell/Research/Mutual Multivariate Independence/Simulation")
} else {
  setwd("~")
}

submit_cran("MDMeasure")
