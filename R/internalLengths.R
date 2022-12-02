#' @title internalLengths
#'
#' @description Provides an estimate for the net growth rate of the clone with
#'     confidence bounds, as well as the sum of internal lengths
#'
#' @param subtree An ape tree subset to include only the clone of interest
#' @param includeStem Boolean indicating whether we should count the stem of the
#'     tree as contributing to the internal lengths summation
#' @param alpha Used for calculation of confidence intervals. 1-alpha confidence
#'    intervals used with default of alpha = 0.05 (95% confidence intervals)
#'
#' @return A named vector including the net growth rate estimate and the sum of
#' internal lengths
#' @examples
#' data(exampleTrees)
#' output_vec <- internalLengths(exampleTrees[[1]])
#' @export
#' @importFrom phangorn "allDescendants"
#' @importFrom ape "is.ultrametric"
internalLengths <- function(subtree, includeStem = F, alpha = 0.05) {

  # Only works for ultrametric trees
  if (! is.ultrametric(subtree)){
    stop("Tree is not ultrametric. internalLengths fn. should only be used with
         ultrametric trees.")
  }

  # Make sure alpha is reasonable
  if(alpha < 0 | alpha > 1){
    stop("alpha must be between 0 and 1")
  }
  if (alpha > 0.25){
    warning(paste0("We calulate 1-alpha confidence intervals. The given confidence
            intervals, with alpha = ", alpha, " correspond to ", 1-alpha, "%
            confidence intervals, which will be very narrow."))
  }

  # Check if tree has stem
  n <- length(subtree$tip.label)
  nodes <- subtree$edge[subtree$edge > n]
  if ( 1 %in% table(nodes)){
    hasStem <- T
    stemNode <- as.numeric(names(which(table(nodes)==1)))
  } else {
    hasStem <- F
  }

  # If includeStem is TRUE, make sure tree has stem
  if(! hasStem & includeStem){
    stop("includeStem is set to TRUE, but tree does not have a stem!")
  }

  # Get list of descendants from each internal node
  descendant_list <- allDescendants(subtree)
  descendant_df <- data.frame("Node" = (length(subtree$tip.label)+2):max(subtree$edge), "Parent" = NA,
                              "Edge_length" = NA, "n_cells" = NA)

  # Find parent, edge length, and number of descendant cells for each internal node
  for (k in descendant_df$Node) {
    descendant_df$n_cells[descendant_df$Node == k] <- sum(descendant_list[[k]] <= length(subtree$tip.label))
    descendant_df$Edge_length[descendant_df$Node == k] <- subtree$edge.length[which(subtree$edge[,2] == k )]
    descendant_df$Parent[descendant_df$Node == k] <- subtree$edge[which(subtree$edge[,2] == k),1]
  }

  # If include Stem is FALSE but tree has a stem, remove stem from calculation
  if ( !includeStem & hasStem){
    descendant_df <- descendant_df[!descendant_df$Parent == stemNode,]
  }

  # The sum of edge lengths in descendant_df is equal to the total internal lengths
  IL <- sum(descendant_df$Edge_length)

  # Calculate growth rate and confidence intervals
  growthRate <- n/IL
  growthRate_lb <- growthRate*(1+qnorm(alpha/2)/sqrt(n))
  growthRate_ub <- growthRate*(1-qnorm(alpha/2)/sqrt(n))

  #return(c(growthRate_lb, growthRate, growthRate_ub, IL))
  return(data.frame("lowerBound" = growthRate_lb, "estimate" = growthRate,
                    "upperBound" = growthRate_ub, "sumInternalLengths" = IL,
                    "n" = n, "alpha" = alpha, "hasStem" = hasStem,
                    "includeStem" = includeStem))
}
