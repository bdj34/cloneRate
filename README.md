
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cloneRate

<!-- badges: start -->
<!-- badges: end -->

The goal of cloneRate is to provide easily accessible methods for
estimating the growth rate of clones using time-based phylogenetic trees
as input. This package provides the internal lengths and maximum
likelihood methods from our recent preprint [Estimating single cell
clonal dynamics in human blood using coalescent
theory](https://www.biorxiv.org)

## Installation

You can install the development version of cloneRate from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bdj34/cloneRate")
```

## Example

### Simulated data

This is a basic example which shows you how to simulate a tree and
calculate the growth rate using two different methods:

``` r
library(cloneRate)

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
#> [1] "Max. likelihood estimate = 0.578961468328489"

# Estimate the growth rate r=a-b=0.5 using internal lengths
intLengths.df <- internalLengths(tree)
print(paste0("Internal lengths estimate = ", intLengths.df$estimate))
#> [1] "Internal lengths estimate = 0.563941151544989"
```

In our [paper](https://www.biorxiv.org), we use simulated trees to test
our growth rate estimates.

### Real data

You can apply these methods to real data as well. Here, we apply our
maximum likelihood and lengths estimates to a real data clone from
[Williams et al.](https://pubmed.ncbi.nlm.nih.gov/35058638/)

``` r
library(cloneRate)

# Load and plot the data
PD9478 <- cloneRate::realCloneData[["PD9478_1"]]
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
print(maxLikelihood(cloneTree)$estimate)
#> [1] 0.5156747
print(internalLengths(cloneTree)$estimate)
#> [1] 0.6010937
```

Our package comes with 42 clones annotated from four distinct
publications, which are the ones we use in our analysis. Note that there
are three clones profiled at two different timepoints, meaning there are
39 unique clones. The papers which generate this data are:

[Williams et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35058638/)

[Mitchell et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35650442/)

[Fabre et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35650444/)

[Van Egeren et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33621486/)
