library(devtools)

# check the source package and record the time
if (.Platform$OS.type == "windows") {
  setwd("C:/Academia/Cornell/Research/Mutual Multivariate Independence/Simulation")
} else {
  setwd("~")
}

system.time(devtools::build_win("MDMeasure", version = "R-release"))
system.time(devtools::build_win("MDMeasure", version = "R-devel"))