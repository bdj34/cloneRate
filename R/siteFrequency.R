#' Get site frequency spectrum of a tree
#'
#' @description `siteFrequency()` calculates the site frequency in
#'  units of time or mutations, as well as a normalized frequency.
#'
#' @param tree An ultrametric or mutation-based tree subset to include only the
#'  clone of interest. Alternatively, a list with several such trees.
#' @param includeStem Boolean indicating whether we should count the stem of the
#'  tree as contributing to the site frequency distribution. Default is FALSE.
#'
#' @returns A data.frame with three columns: the number of descendant cells, site
#'  frequency in units of time or mutations, and normalized site frequency. If a
#'  list of trees is input, output will be a list of such data.frames.
#' @seealso [cloneRate::internalLengths()] and [cloneRate::sharedMuts()] which
#'  use the sum of edge lengths ancestral to between 2 and n-1 tips to calculate
#'  a growth rate.
#' @export
#' @examples
#' # Get site frequency of a single tree
#' example.df <- siteFrequency(exampleUltraTrees[[1]])
#'
#' # Get site frequency of a list of trees
#' example.list <- siteFrequency(exampleMutTrees)
#'
siteFrequency <- function(tree, includeStem = FALSE) {
  # If we have a list of phylo objects instead of a single phylo objects, call recursively
  if (inherits(tree, "list") & !inherits(tree, "phylo")) {
    # Call function recursively on all trees in list
    return.list <- lapply(tree, siteFrequency)

    # Keep naming if they had names
    if (!is.null(names(tree))) {
      names(return.list) <- paste0(names(tree), "_siteFrequency")
    }
    return(return.list)
  }

  # Make sure it's a phylo object
  if (!inherits(tree, "phylo") | is.null(tree$edge.length)) {
    stop("Input tree must be of class 'phylo' (or a list of objects of class
         'phylo') with non-null 'edge.length' vector.")
  }

  if (includeStem) {
    message("You have set includeStem = TRUE. Note that we do not include the stem
            as part of the site frequency calculation in our work (Johnson et
            al. 2022), due to the fact that we don't know when clone initiation
            actually occurs.")
  }

  # Check if tree has stem
  n <- length(tree$tip.label)
  nodes <- tree$edge[tree$edge > n]
  if (1 %in% table(nodes)) {
    hasStem <- TRUE
    stemNode <- as.numeric(names(which(table(nodes) == 1)))
  } else {
    hasStem <- FALSE
  }

  # If includeStem is TRUE, make sure tree has stem
  if (!hasStem & includeStem) {
    stop("includeStem is set to TRUE, but tree does not have a stem!")
  }

  # Get list of descendants from each internal node
  descendant_df <- data.frame(
    "Node" = (length(tree$tip.label) + 2):max(tree$edge), "Parent" = NA,
    "Edge_length" = NA, "n_descendants" = NA
  )


  # Find parent and edge length preceding each internal node
  numDescendants <- getNumberTips(tree)
  for (k in descendant_df$Node) {
    descendant_df$n_descendants[descendant_df$Node == k] <-
      numDescendants$n_tips[numDescendants$Node == k]
    descendant_df$Edge_length[descendant_df$Node == k] <- tree$edge.length[which(tree$edge[, 2] == k)]
    descendant_df$Parent[descendant_df$Node == k] <- tree$edge[which(tree$edge[, 2] == k), 1]
  }

  # If include Stem is FALSE but tree has a stem, remove stem from calculation
  if (!includeStem & hasStem) {
    descendant_df <- descendant_df[!descendant_df$Parent == stemNode, ]
  }

  # Finally, sum up and record the edge lengths which are shared by the same number of cells
  site_freq <- data.frame(
    "n_descendants" = c(2:max(descendant_df$n_descendants)),
    "freq" = NA
  )
  for (k in site_freq$n_descendants) {
    site_freq$freq[site_freq$n_descendants == k] <- sum(descendant_df$Edge_length[descendant_df$n_descendants == k])
  }

  # Add a column for normalized frequency
  site_freq$normalizedFreq <- site_freq$freq / sum(site_freq$freq)

  return(site_freq)
}



#' Get the number of tips that are descendant from an internal node
#'
#' @param tree of class 'phylo'
#'
#' @returns data.frame with one column containing node numbers and the other
#'   containing the number of tips descending from the given node
#'
#' @noRd
getNumberTips <- function(tree) {
  # Define the data.frame matching node to number of descendant tips
  tipNumbers <- data.frame("Node" = c(1:max(tree$edge)), "n_tips" = rep(NA, max(tree$edge)))

  # Get number of tips
  nTips <- ape::Ntip(tree)

  # Set the tips equal to 1
  tipNumbers$n_tips[tipNumbers$Node <= nTips] <- 1

  # Set unknown nodes vector
  unknownNodes <- seq(nTips + 1, max(tree$edge))

  # Fill in the unknown nodes
  while (length(unknownNodes) > 0) {
    for (i in unknownNodes) {
      # Get the two direct descendants of the node
      descendantNodes <- tree$edge[tree$edge[, 1] == i, 2]

      # If we know the number of tips of both direct descendants
      if (!anyNA(tipNumbers$n_tips[tipNumbers$Node %in% descendantNodes])) {
        # sum to find tips of current node
        tipNumbers$n_tips[tipNumbers$Node == i] <-
          sum(tipNumbers$n_tips[tipNumbers$Node %in% descendantNodes])
        # Remove current nodes from unknown nodes
        unknownNodes <- unknownNodes[unknownNodes != i]
      }
    }
  }

  return(tipNumbers)
}
