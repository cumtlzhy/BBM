
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BBM

<!-- badges: start -->

<!-- badges: end -->

The goal of BBM is to mine the co-methylation patterns hidden in the
MeRIP-Seq data. BBM can be applied directly to the IP and input data
matrices for m6A co-methylation mining, avoiding the problem of
introducing noise or information distortion of traditional methods,
which require the IP and input data matrices to be first approximated
into methylation levels. It can also be used for pattern mining of data
that follows a beta-binomial mixture distribution.If you use BBM for
other types of data mining, the input matrices a and b must be numeric
matrices.

## Installation

You can install the released version of BBM from
[github](https://github.com/cumtlzhy/BBM.git) with:

``` r
library(usethis)
library(devtools)
install_github("cumtlzhy/BBM")
```

## Example

This is a basic example which shows you how to use BBM:

``` r

## basic example code
rm(list = ls(all = TRUE))
library(BBM)
data("ip")
data("input")
number_bicluster <- 6
bic_row_labels <- c(rep(0,dim(ip)[1]))
patterns_found_by_BBM <- multiple_GSB(a = ip,b = input,iteration = 1000,burn_in = 500,bic_row_labels = bic_row_labels,number_bicluster = number_bicluster)
```
