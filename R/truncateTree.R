#' Truncate a tree at a distance dist from the root. Useful for identifying subclones
#' and especially for embryonic development rate estimation. Imagine taking a tree, drawing a line across it
#' at a time dist, and only evaluating the tree from above that line.
#'
#' @param tree An ultrametric (time-based) or mutation-based tree of class 'phylo' (see ape package for documentation).
#' @param dist Evolutionary distance (in time or mutations) down the tree at which to truncate.
#' Note that the units should be the same as the units for the tree's edge lengths.
#'
#' @return A truncated tree with the same units of edge lengths as the input tree
#' @export
#'
#' @examples
#' tree <- cloneRate::embryonic_mutation_trees[["KX001"]]
#' dist=55
#' truncated_tree <- truncate_tree(tree, dist)
truncate_tree <- function(tree, dist) {

  # If we have a list of phylo objects instead of a single phylo objects, call recursively
  if (inherits(tree, "list") & !inherits(tree, "phylo")) {
    # Call function recursively on all trees in list, then combine results into one data.frame
    return_list <- lapply(tree, truncate_tree, dist=dist)
    return(return_list)
  }

  # Calculate root to all node distances
  all_node_distances <- ape::node.depth.edgelength(tree)

  # Identify tips with distance less than X
  #tips_to_keep <- names(distances[distances <= X])

  # Identify the nodes with distance less than or equal to dist
  nodes_keep <- which(all_node_distances <= dist)
  nodes_remove <- which(all_node_distances > dist)
  nodes_remove_internal <- nodes_remove[nodes_remove > ape::Ntip(tree)]

  # Identify the nodes with parents who are kept
  parents_are_kept <- c()
  for(i in nodes_remove_internal){
    node <- i
    parent <- getImmediateParent(tree, node)

    if(parent %in% nodes_keep){
      parents_are_kept <- c(parents_are_kept, node)
    }
  }

  # Make sure there's no duplication
  stopifnot(parents_are_kept == unique(parents_are_kept))

  dropTips <- c()
  for(i in parents_are_kept){
    dropTips <- c(dropTips, utils::tail(getTipDescendants(tree, i)[[1]], -1))
  }
  dropTips <- unique(dropTips)

  # Drop all but one tip from each internal node that has depth greater than dist and parent less than dist
  new_tree <- ape::drop.tip(tree, tip = dropTips, collapse.singles = T)

  # Scale the outer edge lengths so every tip has a depth of dist
  index_tip <- which(new_tree$edge[,2] <= ape::Ntip(new_tree))
  new_depths <- ape::node.depth.edgelength(new_tree)
  for(row in index_tip){
    tipNumber <- new_tree$edge[row,2]
    tipDepth <- new_depths[tipNumber]
    newEdgeLength <- new_tree$edge.length[row] - (tipDepth - dist)
    new_tree$edge.length[row] <- newEdgeLength
  }

  return(new_tree)
}






#' Get tip descendants of a tree
#'
#' @description Returns the tips descendant from a node in a phylogenetic tree.
#'
#' @noRd
#' @param tree ape 'phylo' object representing the phylogenetic tree.
#' @param node node number that we want descendants of.
#' @keywords internal
#'
#' @returns numeric. The tip nodes descendant of the input node.
#'
getTipDescendants <- function(tree, node) {
  descendants <- numeric(0)  # Initialize an empty vector to hold descendants
  toVisit <- node  # Start with the node of interest
  nTips <- ape::Ntip(tree)  # Get the number of tip nodes

  while(length(toVisit) > 0) {
    current <- toVisit[1]  # Take the first node to visit
    toVisit <- toVisit[-1]  # Remove it from the list

    # Find direct children of the current node
    children <- tree$edge[tree$edge[,1] == current, 2]

    # Separate children into tips and internal nodes
    tipChildren <- children[children <= nTips]  # Tips have IDs <= nTips
    internalChildren <- children[children > nTips]  # Internal nodes have IDs > nTips

    # Add tip children to the list of descendants
    descendants <- c(descendants, tipChildren)

    # Only add internal nodes to the list of nodes to visit for further exploration
    toVisit <- c(toVisit, internalChildren)
  }

  return(descendants)
}


#' Get direct parent of an ape tree
#'
#' @description Get the immediate ancestor (parent) of a given node for a given tree. Return NA if the node is the root.
#'
#' @noRd
#' @param tree ape 'phylo' object representing the phylogenetic tree.
#' @param node node number that we want the parent of.
#' @keywords internal
#'
#' @returns numeric. The parent node of the input node. NA if input node is root.
#'
getImmediateParent <- function(tree, node) {
  # Look for the node in the second column of the edge matrix to find its parent
  parentRow <- which(tree$edge[,2] == node)

  # If the node has a parent, return the parent node number; otherwise, return NA
  if(length(parentRow) > 0) {
    return(tree$edge[parentRow, 1])
  } else {
    return(NA)  # The node does not have a parent, possibly because it's the root or an error
  }
}
