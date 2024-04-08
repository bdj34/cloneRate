#' Real clone data from human blood
#'
#' 42 clones (39 distinct) from 32 individual donors, 13 of whom have a diagnosis
#' of Myeloproliferative Neoplasm
#'
#' @docType data
#'
#' @usage data(realCloneData)
#'
#' @format A \code{list} of containing one \code{list} with the full ultrametric
#' trees from 30 of the 32 individual donors (the two from Van Egeren are not included),
#' and one \code{list} containing the 42 clone trees.In three cases, there are two
#' timepoints from the same clone, and these are separate phylo objects.
#' Each \code{list} contains a tree as a class \code{phylo} object.
#' See ape package documentation for details on class \code{phylo} objects.
#' Names of each \code{phylo} object (tree) in the list matches the naming used
#' in the sources and also includes driver, age, and clone number.
#'
#'
#' @references These datasets were generated and annotated in:
#' [Williams et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35058638/)
#' [Mitchell et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35650442/)
#' [Fabre et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35650444/)
#' [Van Egeren et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33621486/)
#' @examples
#' # Plot full reconstructed tree from donor PD34493
#' ape::plot.phylo(cloneRate::realCloneData[["fullTrees"]][["PD34493"]],
#'   direction = "downwards", show.tip.label = FALSE
#' )
#'
"realCloneData"


#' Embryonic mutation trees from human blood
#'
#' 12 trees from 12 individual donors.
#'
#' @docType data
#'
#' @usage data(embryonic_mutation_trees)
#'
#' @format A \code{list} of containing trees with mutation-based edge lengths.
#' Two of these trees, CB001 and CB002, are from cord blood samples taken at birth.
#' The remaining trees are from adult samples, truncated to include only the first 55 mutations,
#' based on findings from Mitchell et al. in the cord blood samples. With the exception of
#' the cord blood trees, the truncated adult trees are ultrametric, even though the units
#' of their edge lengths are mutations, and not time.
#' Each element of the \code{list} is a tree as a class \code{phylo} object.
#' See ape package documentation for details on class \code{phylo} objects.
#' Names of each \code{phylo} object (tree) in the list matches the naming used
#' in the sources. Source paper is indicated in the metadata dataframe for each tree.
#'
#'
#' @references These datasets were generated and annotated in:
#' [Williams et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35058638/)
#' [Mitchell et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35650442/)
#' [Fabre et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35650444/)
#' @examples
#' # Plot embryonic mutation tree from cord blood sample CB002
#' ape::plot.phylo(cloneRate::embryonic_mutation_trees[["CB002"]],
#'   direction = "downwards", show.tip.label = FALSE
#' )
#'
#' # Calculate the growth rate using the shared mutations method, assuming 55 mutations per 40/52 weeks (avg. embryonic time of 40 weeks)
#' out <- sharedMutations(cloneRate::embryonic_mutation_trees, nu = 55/(40/52), allow.ultrametric = TRUE)
#'
"embryonic_mutation_trees"


#' Embryonic time-based trees from human blood
#'
#' 12 trees from 12 individual donors.
#'
#' @docType data
#'
#' @usage data(embryonic_time_trees)
#'
#' @format A \code{list} of containing trees with mutation-based edge lengths.
#' Two of these trees, CB001 and CB002, are from cord blood samples taken at birth.
#' The remaining trees are from adult samples, truncated to include only the first 55 mutations,
#' and the scaled to convert mutations to years
#' based on findings from Mitchell et al. in the cord blood samples.
#' Each element of the \code{list} is a tree as a class \code{phylo} object.
#' See ape package documentation for details on class \code{phylo} objects.
#' Names of each \code{phylo} object (tree) in the list matches the naming used
#' in the sources. Source paper is indicated in the metadata dataframe for each tree.
#'
#'
#' @references These datasets were generated and annotated in:
#' [Williams et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35058638/)
#' [Mitchell et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35650442/)
#' [Fabre et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35650444/)
#' @examples
#' # Plot embryonic time-based trees from cord blood sample CB002
#' ape::plot.phylo(cloneRate::embryonic_time_trees[["CB002"]],
#'   direction = "downwards", show.tip.label = FALSE
#' )
#'
#' # Calculate the growth rate using the max likelihood method.
#' out <- maxLikelihood(cloneRate::embryonic_time_trees,  allow.ultrametric = TRUE)
#'
"embryonic_time_trees"


#' Longitudinal validation data
#'
#' For three individuals with clonal expansions that can be estimated using our methods,
#' we have longitudinal data to orthogonally validate these estimates, which is included here.
#' Additionally, for 13 clones with a driver gene matching a driver gene in the single cell
#' data, but without a match to a specific clone, we include this longitudinal data as well.
#'
#' @format A \code{data.frame} containing all the information needed
#' \describe{
#'  \item{Sample.ID}{The individual's ID}
#'  \item{Age}{Individual's age at the various sampling times}
#'  \item{VAF}{The variant allele frequency at the various sampling times for the clone of interest}
#'  \item{Gene}{Gene or genes with mutation that identifies the clone}
#'  \item{Protein}{Protein affected by the mutation}
#'  \item{cellType}{The type of cells used for sequencing}
#'  \item{cloneName}{The name we use for the clone to match to single cell data, if applicable.}
#' }
#'
#' @references These datasets were generated and annotated in:
#' [Williams et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35058638/)
#' [Fabre et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35650444/)
#'
#' @keywords longitudinal
#'
#' @examples
#' # Plot longitudinal data from PD9478
#' library(ggplot2)
#' ggplot(longitudinalData[longitudinalData$Sample.ID == "PD9478", ]) +
#'   geom_point(aes(x = Age, y = VAF))
"longitudinalData"
