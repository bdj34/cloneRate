
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

### Simulated data

This is a basic example which shows you how to simulate a tree and
calculate the growth rate using two different methods:

``` r
library(coalRate)

# Generate a sampled tree with 100 tips from a 20 year birth-death process with birth rate a=1 and death rate b=0.5
tree <- simTree(a=1, b=0.5, cloneAge=20, n=100)

# Plot the tree (see ggtree for more advanced plotting)
plot(tree)
```

<img src="man/figures/README-sim example-1.png" width="100%" />

``` r

# Estimate the growth rate r=a-b=0.5 using maximum likelihood
maxLike.df <- maxLikelihood(tree)
print(paste0("Max. likelihood estimate = ", maxLike.df$estimate))
#> [1] "Max. likelihood estimate = 0.50691856187539"

# Estimate the growth rate r=a-b=0.5 using internal lengths
intLengths.df <- internalLengths(tree)
print(paste0("Internal lengths estimate = ", intLengths.df$estimate))
#> [1] "Internal lengths estimate = 0.521575201726908"
```

In our [paper](https://www.biorxiv.org), we use simulated trees to test
our growth rate estimates.

### Real data

You can apply these methods to real data as well. Here, we apply our
maximum likelihood and lengths estimates to a real data clone from
[Williams et al.](https://pubmed.ncbi.nlm.nih.gov/35058638/)

``` r
library(coalRate)

# Load and plot the data
PD9478 <- coalRate::realCloneData[["PD9478_1"]]
plot(PD9478)
```

<img src="man/figures/README-realData example-1.png" width="100%" />

Let’s now use the expansions data.frame to subset the tree using the ape
function keep.tip()

``` r
# Select the subclone identified in the expansions data.frame
cloneTree <- ape::keep.tip(PD9478, PD9478$expansions$Tips[PD9478$expansions$Expansion == "clone1"])
plot(cloneTree)
```

<img src="man/figures/README-get subclone-1.png" width="100%" />

We can now apply our methods to the clone tree

``` r
# Get maximum likelihood and internal lengths estimates
print(maxLikelihood(cloneTree))
#>   lowerBound  estimate upperBound sumInternalLengths sumExternalLengths
#> 1  0.4153672 0.5156747  0.6159822            118.118           1156.869
#>   extIntRatio  n alpha hasStem includeStem runtime_s  method
#> 1    9.794179 71  0.05   FALSE       FALSE     0.012 maxLike
print(internalLengths(cloneTree))
#>   lowerBound  estimate upperBound sumInternalLengths sumExternalLengths
#> 1  0.4612763 0.6010937   0.740911            118.118           1156.869
#>   extIntRatio  n alpha hasStem includeStem runtime_s  method
#> 1    9.794179 71  0.05   FALSE       FALSE     0.004 lengths
```

Our package comes with 42 clones annotated from four distinct
publications, which are the ones we use in our analysis. Note that there
are three clones profiled at two different timepoints, meaning there are
39 unique clones. The papers which generate this data are:

[Williams et al.](https://pubmed.ncbi.nlm.nih.gov/35058638/)

[Mitchell et al.](https://pubmed.ncbi.nlm.nih.gov/35650442/)

[Fabre et al.](https://pubmed.ncbi.nlm.nih.gov/35650444/)

[Van Egeren et al.](https://pubmed.ncbi.nlm.nih.gov/33621486/)
