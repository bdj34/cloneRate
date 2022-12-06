#' Example tree data
#'
#' Set of 100 trees reconstructed from the distribution of a sample of n=100 tips.
#' All trees have a net growth rate of 1 with birth rates between 1 and 2
#' (sampled from a uniform distribution). Death rates are equal to the chosen
#' birth rate minus 1. Tree reconstruction uses the exact distribution of
#' coalescence times described in "The coalescent of a sample from a binary
#' branching process", Lambert A., Theor. Pop. Bio. 2018. Tree construction and
#' formatting uses \code{ape} R package [ape::rcoal()].
#'
#' @docType data
#'
#' @usage data(exampleTrees)
#'
#' @format A \code{list} of objects of class \code{phylo}
#' \describe{
#'  \item{edge}{A matrix of edge connections which reconstruct the tree.}
#'  \item{edge.length}{A numeric vector of the branch lengths of the connections
#'  in \code{edge} matrix.}
#'  \item{tip.label}{A character vector containing the  (arbitrary in this case)
#'   labels for the 100 tips/samples of the tree.}
#'  \item{Nnode}{Integer number of internal nodes of the tree}
#'  See ape package for details on class \code{phylo} objects.
#'  Names of each \code{phylo} object (tree) in the list provides information
#'  about the net growth rate (r), birth rate (a), death rate (b), number of
#'  tips (n), and clone age (T).
#' }
#' @references This data set was created for the coalRate package using
#' coalescent theory approaches described in "The coalescent of a sample from a
#' binary branching process", Lambert A., Theor. Pop. Bio. 2018.
#' @keywords phylogenetics, trees
#' @examples
#' data(exampleTrees) # Load data
#' library(ggtree) # Package for plotting phylogenetic trees by extending ggplot2
#' ggtree(exampleTrees[[1]]) # Plot first of 100 trees
#'
"exampleTrees"
