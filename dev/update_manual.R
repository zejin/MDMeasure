# create the reference manual
if (.Platform$OS.type == "windows") {
  setwd("C:/Academia/Cornell/Research/Mutual Multivariate Independence")
} else {
  setwd("~")
}

system("R CMD Rd2pdf MDMeasure")
system("mv MDMeasure.pdf MDMeasure/dev/doc")
# system("R CMD check rddapp")
