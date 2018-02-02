library(devtools)

# set up a basic .travis.yml config file
if (.Platform$OS.type == "windows") {
  setwd("C:/Academia/Cornell/Research/Mutual Multivariate Independence")
} else {
  setwd("~")
}

devtools::use_revdep("MDMeasure")
