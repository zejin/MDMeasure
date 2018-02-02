library(devtools)

# resubmit the source package to cran 
if (.Platform$OS.type == "windows") {
  setwd("C:/Academia/Cornell/Research/Mutual Multivariate Independence")
} else {
  setwd("~")
}

submit_cran("MDMeasure")
