
<!-- README.md is generated from README.Rmd. Please edit that file -->

# cloneRate

<!-- badges: start -->
<!-- badges: end -->

The goal of cloneRate is to provide easily accessible methods for
estimating the growth rate of clones. The input should either be an
ultrametric phylogenetic tree with edge lengths corresponding to time,
or a non-ultrametric phylogenetic tree with edge lengths corresponding
to mutation counts. This package provides the internal lengths and
maximum likelihood methods for ultrametric trees and the shared
mutations method for mutation-based trees, all of which are from our
recent preprint [Estimating single cell clonal dynamics in human blood
using coalescent theory](https://www.biorxiv.org)

## Installation

You can install the development version of cloneRate from
[GitHub](https://github.com/) with:

``` r
# Install devtools if you don't have it already
install.packages(setdiff("devtools", rownames(installed.packages()))) 
devtools::install_github("bdj34/cloneRate")
```

For this tutorial, we will also use ggplot2 and ggtree, which can be
installed from CRAN if you don’t have them.

``` r
install.packages(setdiff(c("ggplot2", "ggtree"), rownames(installed.packages())))
```

## Example

### Simulated data

This is a basic example which shows you how to simulate a tree and plot
it.

``` r
library(cloneRate, quietly = T)
library(ggtree, quietly = T)
#> ggtree v3.6.2 For help: https://yulab-smu.top/treedata-book/
#> 
#> If you use the ggtree package suite in published research, please cite
#> the appropriate paper(s):
#> 
#> Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam.
#> ggtree: an R package for visualization and annotation of phylogenetic
#> trees with their covariates and other associated data. Methods in
#> Ecology and Evolution. 2017, 8(1):28-36. doi:10.1111/2041-210X.12628
#> 
#> Shuangbin Xu, Lin Li, Xiao Luo, Meijun Chen, Wenli Tang, Li Zhan, Zehan
#> Dai, Tommy T. Lam, Yi Guan, Guangchuang Yu. Ggtree: A serialized data
#> object for visualization of a phylogenetic tree and annotation data.
#> iMeta 2022, 4(1):e56. doi:10.1002/imt2.56
#> 
#> LG Wang, TTY Lam, S Xu, Z Dai, L Zhou, T Feng, P Guo, CW Dunn, BR
#> Jones, T Bradley, H Zhu, Y Guan, Y Jiang, G Yu. treeio: an R package
#> for phylogenetic tree input and output with richly annotated and
#> associated data. Molecular Biology and Evolution. 2020, 37(2):599-603.
#> doi: 10.1093/molbev/msz240

# Generate a sampled tree with 100 tips from a 20 year birth-death process with birth rate a=1 and death rate b=0.5
tree <- simTree(a=1, b=0.5, cloneAge=30, n=100)

# Plot the tree (see ggtree for more advanced plotting)
ggtree(tree) + layout_dendrogram()
```

<img src="man/figures/README-sim plot-1.png" width="100%" />

Then, we can use this tree as input to our methods for growth rate
estimation:

``` r
# Estimate the growth rate r=a-b=0.5 using maximum likelihood
maxLike.df <- maxLikelihood(tree)
print(paste0("Max. likelihood estimate = ", round(maxLike.df$estimate, 3)))
#> [1] "Max. likelihood estimate = 0.516"

# Estimate the growth rate r=a-b=0.5 using internal lengths
intLengths.df <- internalLengths(tree)
print(paste0("Internal lengths estimate = ", round(intLengths.df$estimate, 3)))
#> [1] "Internal lengths estimate = 0.527"
```

In our [paper](https://www.biorxiv.org), we use simulated trees to test
our growth rate estimates. As an example, let’s load some simulated data
that comes with our package, exampleUltraTrees has 100 ultrametric trees
and exampleMutTrees has 100 mutation-based trees. In the “params”
data.frame we will find the ground truth growth rate, which in this case
is 1. Then, let’s apply our methods to all of these trees.

``` r

# Here we are applying our two methods to all of the ultrametric trees, and then
# using rbind to make one data.frame for each method containing 100 estimates
resultsUltraMaxLike <- do.call(rbind, lapply(exampleUltraTrees, maxLikelihood))
resultsUltraLengths <- do.call(rbind, lapply(exampleUltraTrees, internalLengths))

# Now, let's do the same thing for the shared mutation trees, with the method 
# now being our shared mutations function (analogous to internal Lengths, but with mutations)
resultsMutsShared <- do.call(rbind, lapply(exampleMutTrees, sharedMuts))
```

Now that we have 100 estimates on 100 different trees from 3 different
methods, let’s plot the distributions

``` r
library(ggplot2, quietly = T)

# Combine all into one df for plotting. Note: we have to subset columns because 
# outputs are slightly different between ultrametric and shared mutation trees
common_cols <- intersect(colnames(resultsUltraLengths), colnames(resultsMutsShared))
results.df <- do.call(rbind, list(resultsUltraLengths[,common_cols], 
                  resultsUltraMaxLike[,common_cols], resultsMutsShared[,common_cols]))

# Plot, adding a vertical line at r=1 because that's the true growth rate
ggplot(results.df) + geom_density(aes(x = estimate, color = method)) +
  geom_vline(xintercept = exampleUltraTrees[[1]]$params$r) + theme_bw()
```

<img src="man/figures/README-plotExample-1.png" width="100%" />

Finally, let’s compute the mean, standard deviation, and root mean
square error (RMSE) of the estimates. We expect maximum likelihhod to
perform the best by RMSE, but 100 is a small sample size so anything
could happen…

``` r

mean <- unlist(lapply(split(results.df, results.df$method), 
                      function(x){mean(x$estimate)}))

sd <- unlist(lapply(split(results.df, results.df$method), 
                    function(x){sd(x$estimate)}))

rmse <- unlist(lapply(split(results.df, results.df$method), 
                      function(x){sqrt(sum((x$estimate-x$r)^2)/length(x$estimate))}))

print(mean)
#>    lengths    maxLike sharedMuts 
#>   1.005448   1.005819   1.033520
print(sd)
#>    lengths    maxLike sharedMuts 
#>  0.1151122  0.0887935  0.1083563
print(rmse)
#>    lengths    maxLike sharedMuts 
#>  1.0053067  0.9950473  1.0319770
```

As expected, maximum likelihood performs the best. While the shared
mutations and lengths (also referred to as internalLengths) methods are
based on the same asymptotic result, we’d expect the shared mutation
method to perform slightly worse due to the effect of the randomness in
poissonian mutations. Now, let’s move on to some real data.

### Real data

You can apply these methods to real data as well. Here, we apply our
maximum likelihood and lengths estimates to a real data clone from
[Williams et al.](https://pubmed.ncbi.nlm.nih.gov/35058638/)

``` r
library(cloneRate, quietly = T)
library(ggtree, quietly = T)

# Load and plot the full tree from this individual
PD9478 <- cloneRate::realCloneData$fullTrees$PD9478_1
ggtree(PD9478) + layout_dendrogram() + geom_hilight(node = 85, fill ="red", alpha = .1)
```

<img src="man/figures/README-realData example-1.png" width="100%" />

Let’s now look at the specific clone represented by node 85 (shaded in
red), which we know corresponds to a clone with a JAK2 and DNMT3A
mutation, based on the annotation in [Williams et al. Fig.
3](https://www.nature.com/articles/s41586-021-04312-6)

``` r
# Load the clone (with its overly detailed name)
PD9478_subClone <- cloneRate::realCloneData$cloneTrees[["PD9478_1_clone1_JAK2:p.F537_K539delinsL_AND_DNMT3A:p.Y908*_age68.75_williams"]]

# Plot the clone tree
ggtree(PD9478_subClone)+ layout_dendrogram()
```

<img src="man/figures/README-get subclone-1.png" width="100%" />

We can now apply our methods to the clone tree. We see that this tree is
ultrametric, so we should apply our internalLengths() and
maxLikelihood() functions.

``` r
# Get maximum likelihood and internal lengths estimates
print(maxLikelihood(PD9478_subClone)$estimate)
#> [1] 0.5156747
print(internalLengths(PD9478_subClone)$estimate)
#> [1] 0.6010937
```

Unfortunately, we don’t have a ground truth to compare to when working
with real data. However, this individual clone has longitudinal
sequencing data, which we can use to estimate an orthogonal growth rate.
This will be added shortly.

Our package comes with 42 clones annotated from four distinct
publications, which are the ones we use in our analysis. Note that there
are three clones profiled at two different timepoints, meaning there are
39 unique clones. More vignettes and code to reproduce our full analysis
from our recent [preprint](biorxiv.org) are coming soon! The papers
which generate this data are:

[Williams et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35058638/)
[Mitchell et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35650442/) [Fabre
et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35650444/) [Van Egeren et
al. 2021](https://pubmed.ncbi.nlm.nih.gov/33621486/)
