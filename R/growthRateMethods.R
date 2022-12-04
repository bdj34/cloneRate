#' @title internalLengths
#'
#' @description Provides an estimate for the net growth rate of the clone with
#'     confidence bounds, using the internal lengths method.
#'
#' @param subtree An ape tree subset to include only the clone of interest
#' @param includeStem Boolean indicating whether we should count the stem of the
#'     tree as contributing to the internal lengths summation
#' @param alpha Used for calculation of confidence intervals. 1-alpha confidence
#'    intervals used with default of alpha = 0.05 (95% confidence intervals)
#'
#' @return A dataframe including the net growth rate estimate, the sum of
#' internal lengths and other important details (runtime, n, etc.)
#' @examples
#' data(exampleTrees)
#' df <- internalLengths(exampleTrees[[1]])
#' @export
#' @importFrom phangorn "allDescendants"
#' @importFrom ape "is.ultrametric"
internalLengths <- function(subtree, includeStem = F, alpha = 0.05) {
  ptm <- proc.time()

  # Perform basic checks on the input tree
  inputCheck(subtree, alpha)

  # Check if tree has stem
  n <- length(subtree$tip.label)
  nodes <- subtree$edge[subtree$edge > n]
  if (1 %in% table(nodes)) {
    hasStem <- T
    stemNode <- as.numeric(names(which(table(nodes) == 1)))
  } else {
    hasStem <- F
  }

  # If includeStem is TRUE, make sure tree has stem
  if (!hasStem & includeStem) {
    stop("includeStem is set to TRUE, but tree does not have a stem!")
  }

  # Get list of descendants from each internal node
  descendant_list <- allDescendants(subtree)
  descendant_df <- data.frame(
    "Node" = (length(subtree$tip.label) + 2):max(subtree$edge), "Parent" = NA,
    "Edge_length" = NA, "n_cells" = NA
  )

  # Find parent, edge length, and number of descendant cells for each internal node
  for (k in descendant_df$Node) {
    descendant_df$n_cells[descendant_df$Node == k] <- sum(descendant_list[[k]] <= length(subtree$tip.label))
    descendant_df$Edge_length[descendant_df$Node == k] <- subtree$edge.length[which(subtree$edge[, 2] == k)]
    descendant_df$Parent[descendant_df$Node == k] <- subtree$edge[which(subtree$edge[, 2] == k), 1]
  }

  # If include Stem is FALSE but tree has a stem, remove stem from calculation
  if (!includeStem & hasStem) {
    descendant_df <- descendant_df[!descendant_df$Parent == stemNode, ]
  }

  # The sum of edge lengths in descendant_df is equal to the total internal lengths
  intLen <- sum(descendant_df$Edge_length)

  # Calculate growth rate and confidence intervals
  growthRate <- n / intLen
  growthRate_lb <- growthRate * (1 + qnorm(alpha / 2) / sqrt(n))
  growthRate_ub <- growthRate * (1 - qnorm(alpha / 2) / sqrt(n))

  # Calculate total external lengths
  extLen <- sum(subtree$edge.length[subtree$edge[, 2] %in% c(1:length(subtree$tip.label))])

  # Check ratio of external to internal lengths
  if (extLen / intLen <= 3) {
    warning("External to internal lengths ratio is less than or equal to 3,
            which means internal lengths method may not be applicable.")
  }

  # Get runtime (including all tests)
  runtime <- proc.time() - ptm

  # return(c(growthRate_lb, growthRate, growthRate_ub, intLen))
  return(data.frame(
    "lowerBound" = growthRate_lb, "estimate" = growthRate,
    "upperBound" = growthRate_ub, "sumInternalLengths" = intLen,
    "sumExternalLengths" = extLen,
    "n" = n, "alpha" = alpha, "hasStem" = hasStem,
    "includeStem" = includeStem, "runtime_s" = runtime[["elapsed"]],
    "method" = "lengths"
  ))
}


#' @title moments
#'
#' @description Provides an estimate for the net growth rate of the clone with
#'     confidence bounds using the method of moments.
#'
#' @param subtree An ape tree subset to include only the clone of interest
#' @param alpha Used for calculation of confidence intervals. 1-alpha confidence
#'    intervals used with default of alpha = 0.05 (95% confidence intervals)
#'
#' @return A dataframe including the net growth rate estimate, confidence intervals,
#' and other important details (runtime, n, etc.)
#'
#' @export
#' @importFrom phangorn "allDescendants"
#' @importFrom ape "branching.times"
#' @examples
#' data(exampleTrees)
#' df <- moments(exampleTrees[[1]])
moments <- function(subtree, alpha = 0.05) {
  ptm <- proc.time()

  # Basic check on input formatting and alpha value
  inputCheck(subtree, alpha)

  # Calculate the growth rate
  growthRate <- (pi / sqrt(3)) * 1 / (sd(branching.times(t1)))
  growthRate_lb <- moments_growth_rate * sqrt(1 + 4 * qnorm(alpha / 2) / sqrt(5 * n))
  growthRate_ub <- moments_growth_rate * sqrt(1 - 4 * qnorm(alpha / 2) / sqrt(5 * n))

  # Calculate total internal and external lengths
  extLen <- sum(subtree$edge.length[subtree$edge[, 2] %in% c(1:length(subtree$tip.label))])
  intLen <- internalLengths(subtree, includeStem=F)#$sumInternalLengths


  runtime <- proc.time() - ptm

  return(data.frame(
    "lowerBound" = growthRate_lb, "estimate" = growthRate,
    "upperBound" = growthRate_ub, "sumInternalLengths" = intLen,
    "sumExternalLengths" = extLen,
    "n" = n, "alpha" = alpha, "hasStem" = hasStem,
    "includeStem" = includeStem, "runtime_s" = runtime[["elapsed"]],
    "method" = "moments"
  ))

}


#' @title inputCheck
#'
#' @description Check the validity of inputs to growth rate fns
#'
#' @param subtree An ape tree subset to include only the clone of interest
#' @param alpha Used for calculation of confidence intervals. 1-alpha confidence
#'    intervals used with default of alpha = 0.05 (95% confidence intervals)
#'
#' @return NULL
#'
inputCheck <- function(subtree, alpha){
  # Must be of class phylo
  if (class(subtree) != "phylo") {
    stop("Tree must be of class phylo. Use as.phylo function to convert if the
    formatting is correct. Otherwise, see ape documentation
    https://cran.r-project.org/web/packages/ape/ape.pdf")
  }

  # Only works for ultrametric trees
  if (!is.ultrametric(subtree)) {
    stop("Tree is not ultrametric. intenralLengths, moments, and maxLike fns.
        should only be used with ultrametric trees. For usage with mutation trees,
        use sharedMuts fn.")
  }

  # Make sure alpha is reasonable
  if (alpha < 0 | alpha > 1) {
    stop("alpha must be between 0 and 1")
  }
  if (alpha > 0.25) {
    warning(paste0("We calulate 1-alpha confidence intervals. The given confidence
            intervals, with alpha = ", alpha, " correspond to ", 1 - alpha, "%
            confidence intervals, which will be very narrow."))
  }

  return(NULL)
}


