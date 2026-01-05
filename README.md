
<!-- README.md is generated from README.Rmd. Please edit that file -->

# inferRecom

<!-- badges: start -->

<!-- badges: end -->

inferRecom provides tools for analyzing meiotic recombination events and
runs of homozygosity in family-based genetic data.

## Installation

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("inferRecom")
```

## Usage

To detect recombinations from maternal meioses in a family with 3+
children:

``` r
library(inferRecom)

# Detect crossovers
crossovers <- xoDetect(
  plinkFile = "path/to/data",
  mapFile = "path/to/map.txt",
  familySize = 3,
  parent = "mother"
)
```

## Vignette

A detailed vignette is available with the package to demonstrate the use
of each of the functions.

## Citation

If you use this package, please cite:


    Citation information

## Getting Help

Bug reports and feature requests:
<https://github.com/yourusername/inferRecom/issues>
