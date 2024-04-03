
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
#' tree <- adjusted_mut_tree_list[["KX001"]]
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
    parent <- Ancestors(tree, node, type = "parent")
    # Make sure sibling isn't already included
    #sibling <- Siblings(tree, node)
    #nonTip_siblings <- sibling[sibling > ape::Ntip(tree)]
    if(parent %in% nodes_keep){
      parents_are_kept <- c(parents_are_kept, node)
    }
  }

  # Make sure there's no duplication
  stopifnot(parents_are_kept == unique(parents_are_kept))

  dropTips <- c()
  for(i in parents_are_kept){
    dropTips <- c(dropTips, tail(Descendants(tree, i, type = "tips")[[1]], -1))
  }
  dropTips <- unique(dropTips)

  # Drop all but one tip from each internal node that has depth greater than dist and parent less than dist
  new_tree <- drop.tip(tree, tip = dropTips, collapse.singles = T)

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





#' Get a maximum likelihood estimate of the embryonic rate of expansion from a full
#' tree measured in mutations
#'
#' @param tree An ultrametric (time-based) or mutation-based tree of class 'phylo' (see ape package for documentation).
#' @param dist Evolutionary distance (in time or mutations) down the tree at which to truncate.
#' Note that the units should be the same as the units for the tree's edge lengths.
#'
#' @return A truncated tree with the same units of edge lengths as the input tree
#' @noRd
#' @examples
embryonic_maxLike <- function(tree, n_muts_avg = 55, n_muts_truncate = 100, cloneAge = 40/52) {

  # If we have a list of phylo objects instead of a single phylo objects, call recursively
  if (inherits(tree, "list") & !inherits(tree, "phylo")) {
    # Call function recursively on all trees in list, then combine results into one data.frame
    return.df <- do.call(rbind, lapply(tree, embryonic_maxLike))
    return.df$cloneName_result <- names(tree)
    return(return.df)
  }

  # Calculate root to all node distances
  all_node_distances <- ape::node.depth.edgelength(tree)

  # Identify the nodes with distance greater than dist
  nodes_keep <- which(all_node_distances <= n_muts_truncate)

  # Find the difference between the latest early coalescence and the next coalescence (if it exists)
  if(length(nodes_keep)/length(all_node_distances) > 0.9){
    agedf <- data.frame("age" = cloneAge, "tip.label" = tree$tip.label)

    tree$agedf <- agedf

    tree$edge.length <- round(tree$edge.length)

    ultra_tree <- make.ultrametric.tree(tree)
    return(maxLikelihood(ultra_tree))
    #print("same")
    #return(sharedMuts(tree, nu = n_muts_avg/cloneAge))
  }

  if(length(nodes_keep) > 10){
    coalTree <- coal_to_tree(cloneAge*(n_muts_truncate - all_node_distances[nodes_keep])/n_muts_avg)
    #print(ggtree(coalTree, layout = "dend"))
    return(maxLikelihood(coalTree))
  }else {
    return(NULL)
  }
}



#' Get a sharedMutations estimate of the embryonic rate of expansion from a full
#' tree measured in mutations
#'
#' @param tree An ultrametric (time-based) or mutation-based tree of class 'phylo' (see ape package for documentation).
#' @param dist Evolutionary distance (in time or mutations) down the tree at which to truncate.
#' Note that the units should be the same as the units for the tree's edge lengths.
#'
#' @return A truncated tree with the same units of edge lengths as the input tree
#' @noRd
#'
#' @examples
embryonic_sharedMuts <- function(tree, n_muts_avg = 55, n_muts_truncate = 120, cloneAge = 40/52, alpha = 0.05) {

  # If we have a list of phylo objects instead of a single phylo objects, call recursively
  if (inherits(tree, "list") & !inherits(tree, "phylo")) {
    # Call function recursively on all trees in list, then combine results into one data.frame
    return.df <- do.call(rbind, lapply(tree, embryonic_sharedMuts))
    return.df$cloneName_result <- names(tree)
    return(return.df)
  }

  # Calculate root to all node distances
  all_node_distances <- ape::node.depth.edgelength(tree)

  # Identify the nodes with distance greater than dist
  nodes_keep <- which(all_node_distances <= n_muts_truncate)
  internal_nodes_keep <- nodes_keep[nodes_keep > ape::Ntip(tree)+2]

  if(length(nodes_keep) > 10){
    # Get list of descendants from each internal node
    descendant_df <- data.frame(
      "Node" = internal_nodes_keep, "Parent" = NA,
      "Edge_length" = NA, "n_cells" = NA
    )


    # Find parent, edge length, and number of descendant cells for each internal node
    for (k in descendant_df$Node) {
      descendants <- tree$edge[tree$edge[, 1] == k, 2]
      descendant_df$n_cells[descendant_df$Node == k] <- length(descendants)
      descendant_df$Edge_length[descendant_df$Node == k] <- tree$edge.length[which(tree$edge[, 2] == k)]
      descendant_df$Parent[descendant_df$Node == k] <- tree$edge[which(tree$edge[, 2] == k), 1]
    }

    # If tree has a stem, remove stem from calculation
    #if (hasStem) {
    #  descendant_df <- descendant_df[!descendant_df$Parent == stemNode, ]
    #}

    # The sum of edge lengths in descendant_df is equal to the total shared mutations
    sharedMutations <- sum(descendant_df$Edge_length)

    nu <- n_muts_avg/cloneAge
    n <- length(internal_nodes_keep) + 3
    growthRate <- n * nu / sharedMutations
    print(nu)
    print(n)
    print(sharedMutations)
    lb <- growthRate * (1 + (stats::qnorm(alpha / 2) / sqrt(n)) * (1 + n / sharedMutations))
    ub <- growthRate * (1 - (stats::qnorm(alpha / 2) / sqrt(n)) * (1 + n / sharedMutations))
    # return data.frame
    result.df <- data.frame(
      "lowerBound" = lb, "estimate" = growthRate,
      "upperBound" = ub, "nu" = nu,
      "sharedMutations" = sharedMutations,
      "n" = n, "alpha" = alpha,
      "method" = "sharedMuts"
    )

    return(result.df)
  }else {
    return(NULL)
  }
}

