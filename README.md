
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cloneRate <a href="https://bdj34.github.io/cloneRate/"><img src="man/figures/logo.png" align="right" height="136" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/bdj34/cloneRate/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bdj34/cloneRate/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/bdj34/cloneRate/branch/main/graph/badge.svg)](https://app.codecov.io/gh/bdj34/cloneRate?branch=main)
<!-- badges: end -->

The goal of cloneRate is to provide easily accessible methods for
estimating the growth rate of clones. The input should either be an
ultrametric phylogenetic tree with edge lengths corresponding to time,
or a non-ultrametric phylogenetic tree with edge lengths corresponding
to mutation counts. In the case of mutation-based edge lengths, a
mutation rate estimate should also be provided. This package provides
the internal lengths and maximum likelihood methods for ultrametric
trees and the shared mutations method for mutation-based trees, all of
which are from our recent preprint [Estimating single cell clonal
dynamics in human blood using coalescent
theory](https://www.biorxiv.org)

## Installation

You can install the development version of cloneRate from
[GitHub](https://github.com/) with:

``` r
# Install devtools if you don't have it already
install.packages(setdiff("devtools", rownames(installed.packages())))

# Install 
devtools::install_github("bdj34/cloneRate")
```

For this basic tutorial and our vignettes, we will also use several
other packages, which can all be installed from CRAN. Because these are
listed as packages we suggest, running the following command will
install them along with the vignettes.

``` r
devtools::install_github("bdj34/cloneRate", build_vignettes = TRUE, dependencies = TRUE)
```

Alternatively, you can install them manually:

``` r
install.packages(setdiff(c("ggplot2", "ggtree", "ggsurvfit", "survival", "car"), rownames(installed.packages())))
```

## Example

We’ll walk through simulating a single tree and plotting it, then apply
our growth rate methods.

### Simulate data

We can simulate a sample of size n from a birth-death tree as follows:

``` r
library(cloneRate, quietly = T)
library(ggtree, quietly = T)
library(ggplot2, quietly = T)

# Generate a sampled tree with 100 tips from a 20 year birth-death process with birth rate a=1 and death rate b=0.5
tree <- simUltra(a = 1, b = 0.5, cloneAge = 40, n = 100)
```

Now that we have simulated the tree, let’s plot it:

``` r
# Plot the tree (see ggtree docs for more advanced plotting)
ggtree(tree) + theme_tree2(fgcolor = "blue") + xlab("Time (years)")
```

<img src="man/figures/README-plotTree-1.png" width="100%" />

The `theme_tree2()` function in `ggtree` allows for a time scale to be
added to the tree.

### Estimate growth rate of one tree

We can use this tree as input to our methods for growth rate estimation:

``` r
# Estimate the growth rate r=a-b=0.5 using maximum likelihood
maxLike.df <- maxLikelihood(tree)
print(paste0("Max. likelihood estimate = ", round(maxLike.df$estimate, 3)))
#> [1] "Max. likelihood estimate = 0.452"

# Estimate the growth rate r=a-b=0.5 using internal lengths
intLengths.df <- internalLengths(tree)
print(paste0("Internal lengths estimate = ", round(intLengths.df$estimate, 3)))
#> [1] "Internal lengths estimate = 0.506"
```

### Estimate growth rate of many trees

In our [paper](https://www.biorxiv.org), we use simulated trees to test
our growth rate estimates. As an example, let’s load some simulated data
that comes with our package, exampleUltraTrees has 100 ultrametric
trees. In the “metadata” data.frame we will find the ground truth growth
rate, which in this case is 1. Let’s apply our methods to all 100 trees.

``` r

# Here we are applying two methods to all of the ultrametric trees
resultsUltraMaxLike <- maxLikelihood(exampleUltraTrees)
resultsUltraLengths <- internalLengths(exampleUltraTrees)
```

Notice how the functions `maxLikelihood()` and `internalLengths()` can
take as input either a single tree or a list of trees. Either way, the
output is a `data.frame` containing the results. Now that we have 100
estimates on 100 different trees from 2 different methods, let’s plot
the distributions

``` r
library(ggplot2, quietly = T)

# Combine all into one df for plotting. This works because the columns are the same
resultsCombined <- rbind(resultsUltraMaxLike, resultsUltraLengths)

# Plot, adding a vertical line at r=1 because that's the true growth rate
ggplot(resultsCombined) +
  geom_density(aes(x = estimate, color = method)) +
  geom_vline(xintercept = exampleUltraTrees[[1]]$metadata$r) +
  theme_bw()
```

<img src="man/figures/README-plotExample-1.png" width="100%" />

Finally, let’s compute the root mean square error (RMSE) of the
estimates. We expect maximum likelihhod to perform the best by RMSE, but
100 is a small sample size so anything could happen…

``` r
# Calculate the RMSE
groundTruth <- exampleUltraTrees[[1]]$metadata$r[1]
rmse <- unlist(lapply(
  split(resultsCombined, resultsCombined$method),
  function(x) {
    sqrt(sum((x$estimate - groundTruth)^2) / length(x$estimate))
  }
))

print(rmse)
#>    lengths    maxLike 
#> 0.10049503 0.08869756
```

As expected, maximum likelihood performs the best. Note that this may
change if we regenerate the data. For more details, see our vignettes:

``` r
vignette("cloneRate-dataAnalysis")
vignette("cloneRate-simulate")
```

## References

Our package comes with 42 clones annotated from four distinct
publications, which are the ones we use in our analysis. Note that there
are three clones profiled at two different timepoints, meaning there are
39 unique clones. More vignettes and code to reproduce our full analysis
from our recent [preprint](biorxiv.org) are coming soon! The papers
which generate this data are:

- [Williams et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35058638/)
- [Mitchell et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35650442/)
- [Fabre et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35650444/)
- [Van Egeren et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33621486/)

The mathematical basis for our estimates is detailed in full in [our
paper](https://www.biorxiv.org/).

Simulating the birth-death trees is a direct result of the work of
Amaury Lambert in:

- [Lambert, 2018](https://pubmed.ncbi.nlm.nih.gov/29704514/)
