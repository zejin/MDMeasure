library(devtools)

# generate documents in .Rd from comments
if (.Platform$OS.type == "windows") {
  setwd("C:/Academia/Cornell/Research/Mutual Multivariate Independence")
} else {
  setwd("~")
}

devtools::document("MDMeasure")
cat('useDynLib(MDMeasure, .registration = TRUE)', 
    file = 'MDMeasure/NAMESPACE', append = TRUE)

