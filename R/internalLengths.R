#' @title internalLengths
#'
#' @description Provides an estimate for the net growth rate of the clone as
#'     well as the sum of internal lengths
#'
#' @param subtree An ape tree subset to include only the clone of interest
#' @param withZeroes Boolean indicating whether the subtree includes a zeroes
#'     outgroup as a tip
#'
#' @return A named vector including the net growth rate estimate and the sum of
#' internal lengths
#' @examples
#' data(exampleTrees)
#' output_vec <- internalLengths(exampleTrees[[1]])
#' @export
#' @importFrom phangorn "allDescendants"
#' @importFrom ape "dist.nodes"
internalLengths <- function(subtree, withZeroes = F) {

  # Only works for ultrametric trees which remove original zeros root
  # or ones where we count edge from original zeros root

  #(nu*n)/(M_i)
  # Which becomes, given an ultrametric tree...
  # n/(internal branch times)


  # Use subset tree to get the edge lengths corresponding to the number of descendants
  # First get list of descendents from each internal node
  descendant_list <- allDescendants(subtree)
  descendant_df <- data.frame("Node" = (length(subtree$tip.label)+2):max(subtree$edge), "Parent" = NA,
                              "Edge_length" = NA, "n_cells" = NA)

  # Then find parent and edge length corresponding to each node
  dist_node_mat <- dist.nodes(subtree)
  for (k in descendant_df$Node) {
    descendant_df$n_cells[descendant_df$Node == k] <- sum(descendant_list[[k]] < length(subtree$tip.label)+1)
    descendant_df$Parent[descendant_df$Node == k] <- subtree$edge[which(subtree$edge[,2] == k ),1]
    descendant_df$Edge_length[descendant_df$Node == k] <- dist_node_mat[descendant_df$Parent[descendant_df$Node == k], k]
    # Check edge length to make sure it's the same using dist.nodes or tree itself
    stopifnot(descendant_df$Edge_length[descendant_df$Node == k] == subtree$edge.length[which(subtree$edge[,2] == k )])
  }

  #write.csv(descendant_df, paste0(outDir, "/descendant_info_",
  #                                patient_name, "_age_", round(age,2), "_driver_", gsub(":", "_", driver), ".csv"))

  # Finally, sum up and record the edge lengths which are shared by the same number of cells
  site_freq <- data.frame("n_cells" = unique(descendant_df$n_cells), "freq" = NA)
  for (k in site_freq$n_cells) {
    site_freq$freq[site_freq$n_cells == k] <- sum(descendant_df$Edge_length[descendant_df$n_cells == k])
  }

  # The sum of edge lengths in site_freq dataframe should equal the total internal branch lengths
  site_freq_internal <- sum(site_freq$freq)

  ######## FIX: Get the number of actual tips (exclude zeroes outgroup)
  if (withZeroes) {
    growthRate <- (length(subtree$tip.label) - 1)/site_freq_internal
  } else {
    growthRate <- length(subtree$tip.label)/site_freq_internal
  }

  return(c(growthRate, site_freq_internal))
}
