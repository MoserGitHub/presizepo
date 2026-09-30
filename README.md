
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `presizepo`

`presizepo` is a precision-based (confidence-interval width) sample size
approach for the log odds ratio in a proportional-odds model, for either
a binary or continuouse predictor

## Installation

`presizepo` can be installed via

``` r
install.packages('presizepo', repos = c('https://dcr-unibe-ch.r-universe.dev', 'https://cloud.r-project.org'))
```

or via

``` r
remotes::install_github("dcr-unibe-ch/presizepo")
```

This may require `Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")`
if packages were built under a different R version to the one you are
using.

## Usage

``` r
library(presizepo)
```
