library(devtools)

# install a package from the source package
if (.Platform$OS.type == "windows") {
  setwd("C:/Academia/Cornell/Research/Mutual Multivariate Independence")
} else {
  setwd("~")
}

devtools::install("MDMeasure", build_vignettes = TRUE)
