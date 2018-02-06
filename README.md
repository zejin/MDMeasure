# MDMeasure

[![Travis-CI Build Status](https://travis-ci.org/zejin/MDMeasure.svg?branch=master)](https://travis-ci.org/zejin/MDMeasure.svg?branch=master)

## Overview

**MDMeasure** provides measures of mutual dependence and tests of mutual independence. 

The two main parts are:
- measuring mutual dependence
- testing mutual independence

## Measuring mutual dependence

The mutual dependence measures include:
- asymmetric measure based on distance covariance
- symmectric measure based on distance covariance
- complete measure based on complete V-statistics
- simplified complete measure based on incomplete V-statistics
- asymmetric measure based on complete measure
- simplified asymmetric measure based on simplified complete measure
- symmectric measure based on complete measure
- simplified symmectric measure based on simplified complete measure

## Testing mutual independence

The mutual independence tests based on the mutual dependence measures are implemented as permutation tests.

## Installation

``` r
# Install the released version from CRAN
install.packages("MDMeasure")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("zejin/MDMeasure")
```




