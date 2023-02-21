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
#' @keywords phylogenetics, hematopoiesis
#' @examples
#' # Plot full reconstructed tree from donor PD34493
#' ape::plot.phylo(cloneRate::realCloneData[["fullTrees"]][["PD34493"]],
#'  direction = "downwards", show.tip.label = FALSE)
#'
#' @source <https://pubmed.ncbi.nlm.nih.gov/35058638/>
#' @source <https://pubmed.ncbi.nlm.nih.gov/35650442/>
#' @source <https://pubmed.ncbi.nlm.nih.gov/35650444/>
#' @source <https://pubmed.ncbi.nlm.nih.gov/33621486/>
"realCloneData"







#' Example ultrametric tree data
#'
#' Set of 100 time-based ultrametric trees reconstructed from the distribution
#' of a sample of n=100 tips.
#' All trees have a net growth rate of 1 with birth rates between 1 and 2
#' (sampled from a uniform distribution). Death rates are equal to the chosen
#' birth rate minus 1. Tree reconstruction uses the exact distribution of
#' coalescence times described in "The coalescent of a sample from a binary
#' branching process", Lambert A., Theor. Pop. Bio. 2018. Tree construction and
#' formatting uses \code{ape} R package [ape::rcoal()].
#'
#' @docType data
#'
#' @usage data(exampleUltraTrees)
#'
#' @format A \code{list} of objects of class \code{phylo}
#' \describe{
#'  \item{edge}{A matrix of edge connections which reconstruct the tree.}
#'  \item{edge.length}{A numeric vector of the branch lengths of the connections
#'  in \code{edge} matrix. Units are years.}
#'  \item{tip.label}{A character vector containing the  (arbitrary in this case)
#'   labels for the 100 tips/samples of the tree.}
#'  \item{Nnode}{Integer number of internal nodes of the tree}
#'  \item{params}{\code{data.frame} containing info on the params used to generate
#'  the tree}
#'  See ape package for details on class \code{phylo} objects.
#' }
#' @references This data set was created for the cloneRate package using
#' coalescent theory approaches described in "The coalescent of a sample from a
#' binary branching process", Lambert A., Theor. Pop. Bio. 2018.
#' @keywords phylogenetics, birth-death trees, Coalescent Point Process.
#' @examples
#' # Plot first of 100 trees
#' ape::plot.phylo(cloneRate::exampleUltraTrees[[1]],
#'   direction = "downwards", show.tip.label = FALSE)
"exampleUltraTrees"







#' Example mutation tree data
#'
#' Set of 100 mutation based trees reconstructed from the distribution
#' of a sample of n=100 tips.
#' All trees have a net growth rate of 1 with birth rates between 1 and 2
#' (sampled from a uniform distribution). Death rates are equal to the chosen
#' birth rate minus 1. Tree reconstruction uses the exact distribution of
#' coalescence times described in "The coalescent of a sample from a binary
#' branching process", Lambert A., Theor. Pop. Bio. 2018. Tree construction and
#' formatting uses \code{ape} R package [ape::rcoal()]. We then change the edge
#' lengths from time-based to mutation-based by drawing from a poisson
#' distribution with mean equal to edge length (in units of time) multiplied
#' by the mutation rate, nu, which is drawn from a uniform distribution between
#' 10 and 20 mutations per year.
#'
#' @docType data
#'
#' @usage data(exampleMutTrees)
#'
#' @format A \code{list} of objects of class \code{phylo}
#' \describe{
#'  \item{edge}{A matrix of edge connections which reconstruct the tree.}
#'  \item{edge.length}{A numeric vector of the branch lengths of the connections
#'  in \code{edge} matrix. Units are mutations.}
#'  \item{tip.label}{A character vector containing the  (arbitrary in this case)
#'   labels for the 100 tips/samples of the tree.}
#'  \item{Nnode}{Integer number of internal nodes of the tree}
#'  \item{params}{\code{data.frame} containing info on the params used to generate
#'   the tree}
#'  See ape package for details on class \code{phylo} objects.
#' }
#' @references This data set was created for the cloneRate package using
#' coalescent theory approaches described in "The coalescent of a sample from a
#' binary branching process", Lambert A., Theor. Pop. Bio. 2018.
#' @keywords phylogenetics, birth-death trees, Coalescent Point Process
#' @examples
#' # Plot first of 100 trees
#' ape::plot.phylo(cloneRate::exampleMutTrees[[1]],
#'  direction = "downwards", show.tip.label = FALSE)
#'
"exampleMutTrees"





#' Longitudinal validation data
#'
#' For three individuals with clonal expansions that can be estimated using our methods,
#' we have longitudinal data to orthogonally validate these estimates, which is provided here.
#'
#' @format A \code{data.frame} containing all the information needed
#' \describe{
#'  \item{Sample.ID}{The individual's ID}
#'  \item{Age}{Individual's age at the various sampling times}
#'  \item{VAF}{The variant allele frequency at the various sampling times for the clone of interest}
#'  \item{Gene}{Gene or genes with mutation that identifies the clone}
#'  \item{Protein}{Protein affected by the mutation}
#'  \item{cellType}{The type of cells used for sequencing}
#'  \item{cloneName}{The name we use for the clone to match to single cell data}
#' }
#'
#' @references These datasets were generated and annotated in:
#' [Williams et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35058638/)
#' [Fabre et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35650444/)
#'
#' @keywords longitudinal
#'
#' @examples
#' library(ggplot2)
#' # Plot longitudinal data from PD9478
#' ggplot(longitudinalData[longitudinalData$Sample.ID == "PD9478", ]) +
#'   geom_point(aes(x = Age, y = VAF))
"longitudinalData"
