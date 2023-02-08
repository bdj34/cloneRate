
<!-- README.md is generated from README.Rmd. Please edit that file -->

# coalRate

<!-- badges: start -->
<!-- badges: end -->

The goal of coalRate is to provide easily accessible methods for
estimating the growth rate of clones using time-based phylogenetic trees
as input. This package provides the internal lengths and maximum
likelihood methods from our recent preprint [Estimating single cell
clonal dynamics in human blood using coalescent
theory](https://www.biorxiv.org)

## Installation

You can install the development version of coalRate from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bdj34/coalRate")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(coalRate)

# Generate a sampled tree with 100 tips from a 20 year birth-death process with birth rate a=1 and death rate b=0.5
tree <- simTree(a=1, b=0.5, cloneAge=20, n=100)

# Estimate the growth rate r=a-b=0.5 using maximum likelihood
result.df <- maxLikelihood(tree)
print(result.df$estimate)
#> [1] 0.5488015
```

You can apply these methods to real data as well. Here, we apply our
maximum likelihood and lengths estimates to a real data clone from
[Williams et al.](https://pubmed.ncbi.nlm.nih.gov/35058638/)

``` r
library(coalRate)
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
