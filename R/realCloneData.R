#' Real clone data from human blood
#'
#' 42 clones (39 distinct) from 32 individual donors, 13 of whom have a diagnosis
#' of Myeloproliferative Neoplasm
#'
#' @docType data
#'
#' @usage data(realCloneData)
#'
#' @format A \code{list} of objects of class \code{phylo} with an included
#' \code{data.frame} "expansions" which has annotated clonal expansions
#' \describe{
#'  \item{edge}{A matrix of edge connections which reconstruct the tree.}
#'  \item{edge.length}{A numeric vector of the branch lengths of the connections
#'  in \code{edge} matrix.}
#'  \item{tip.label}{A character vector containing the  (arbitrary in this case)
#'   labels for the 100 tips/samples of the tree.}
#'  \item{Nnode}{Integer number of internal nodes of the tree}
#'  \item{age}{Numeric age of the person at time of sampling}
#'  \item{has_expanded_clades}{Logical indicating whether or not the tree has
#'  somatic clonal expansions of at least 10 sampled cells}
#'  \item{expansions}{data.frame identifying the tips that are to be included in
#'  analysis of expanded clones/clades.}
#'  See ape package for details on class \code{phylo} objects.
#'  Names of each \code{phylo} object (tree) in the list matches the naming used
#'  in the sources.
#' }
#' @references These datasets were generated and annotated in:
#' [Williams et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35058638/)
#' [Mitchell et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35650442/)
#' [Fabre et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35650444/)
#' [Van Egeren et al. 2021](https://pubmed.ncbi.nlm.nih.gov/33621486/)
#' @keywords phylogenetics, hematopoiesis
#' @examples
#' plot(coalRate::realCloneData[["PD34493"]]) # Plot full reconstructed tree from donor PD34493
#'
#'
#' @source <https://pubmed.ncbi.nlm.nih.gov/35058638/>
#' @source <https://pubmed.ncbi.nlm.nih.gov/35650442/>
#' @source <https://pubmed.ncbi.nlm.nih.gov/35650444/>
#' @source <https://pubmed.ncbi.nlm.nih.gov/33621486/>
"realCloneData"
