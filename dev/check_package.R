library(devtools)

# check the source package and record the time
if (.Platform$OS.type == "windows") {
  setwd("C:/Academia/Cornell/Research/Mutual Multivariate Independence")
} else {
  setwd("~")
}

system.time(devtools::check("MDMeasure", document = FALSE, manual = TRUE))
