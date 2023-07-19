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
#' @examples
#' # Plot first of 100 trees
#' ape::plot.phylo(cloneRate::exampleUltraTrees[[1]],
#'   direction = "downwards", show.tip.label = FALSE
#' )
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
#' @examples
#' # Plot first of 100 trees
#' ape::plot.phylo(cloneRate::exampleMutTrees[[1]],
#'   direction = "downwards", show.tip.label = FALSE
#' )
#'
"exampleMutTrees"
