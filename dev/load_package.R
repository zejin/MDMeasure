library(devtools)

# load a package into memory from the source package
if (.Platform$OS.type == "windows") {
  setwd("C:/Academia/Cornell/Research/Mutual Multivariate Independence")
} else {
  setwd("~")
}

devtools::load_all("MDMeasure")

