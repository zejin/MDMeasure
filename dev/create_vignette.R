library(devtools)

# creare a vignette
if (.Platform$OS.type == "windows") {
  setwd("C:/Academia/Cornell/Research/Mutual Multivariate Independence")
} else {
  setwd("~")
}

devtools::use_vignette("MDMeasure", pkg = "MDMeasure")
