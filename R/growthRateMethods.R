#' Growth rate estimate using the sum of internal lengths
#'
#' @description `internalLengths()` provides an estimate for the net growth rate of the clone with confidence bounds, using the internal lengths method.
#'
#' @param subtree An ultrametric tree subset to include only the clone of
#' interest. Alternatively, a list with several such trees.
#' @param includeStem Boolean indicating whether we should count the stem of the tree as contributing to the internal lengths summation
#' @param alpha Used for calculation of confidence intervals. 1-alpha confidence intervals used with default of alpha = 0.05 (95 percent confidence intervals)
#'
#' @returns A dataframe including the net growth rate estimate, the sum of internal lengths and other important details (runtime, n, etc.)
#' @seealso [cloneRate::maxLikelihood()], [cloneRate::sharedMuts()]
#' @export
#' @examples
#' internalLengths(cloneRate::exampleUltraTrees[[1]])
#'
internalLengths <- function(subtree, includeStem = F, alpha = 0.05) {
  ptm <- proc.time()

  # If we have a list of phylo objects instead of a single phylo objects, call recursively
  if (inherits(subtree, "list") & !inherits(subtree, "phylo")) {
    # Call function recursively on all trees in list, then combine results into one data.frame
    return.df <- do.call(rbind, lapply(subtree, internalLengths))
    return.df$names <- names(subtree)
    return(return.df)
  }

  # Perform basic checks on the input tree
  inputCheck(subtree, alpha)

  if (includeStem) {
    message("You have set includeStem = T. Note that we do not include the stem
            as part of the internal lengths calculation in our work (Johnson et al. 2022)")
  }

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
  descendant_df <- data.frame(
    "Node" = (length(subtree$tip.label) + 2):max(subtree$edge), "Parent" = NA,
    "Edge_length" = NA
  )


  # Find parent and edge length preceding each internal node
  for (k in descendant_df$Node) {
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
  growthRate_lb <- growthRate * (1 + stats::qnorm(alpha / 2) / sqrt(n))
  growthRate_ub <- growthRate * (1 - stats::qnorm(alpha / 2) / sqrt(n))

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
    "sumExternalLengths" = extLen, extIntRatio = extLen / intLen,
    "n" = n, "alpha" = alpha, "hasStem" = hasStem,
    "includeStem" = includeStem, "runtime_s" = runtime[["elapsed"]],
    "method" = "lengths"
  ))
}



#' Growth rate estimate using the sum of shared mutations assuming a mutation tree
#'
#' @description `sharedMuts()` provides an estimate for the net growth rate of the clone with confidence bounds, using the shared mutations method.
#'
#' @param subtree A non-ultrametric ape tree subset to include only the clone of interest
#' @param nu The mutation rate
#' @param includeStem Boolean indicating whether we should count the stem of the tree as contributing to the internal lengths summation
#' @param alpha Used for calculation of confidence intervals. 1-alpha confidence intervals used with default of alpha = 0.05 (95 percent confidence intervals)
#'
#' @returns A dataframe including the net growth rate estimate, the sum of internal lengths and other important details (runtime, n, etc.)
#' @seealso [cloneRate::internalLengths()]
#' @export
#' @examples
#' sharedMuts(cloneRate::exampleMutTrees[[1]])
#'
sharedMuts <- function(subtree, nu = NULL, includeStem = F, alpha = 0.05) {
  ptm <- proc.time()

  # If we have a list of phylo objects instead of a single phylo objects, call recursively
  if (inherits(subtree, "list") & !inherits(subtree, "phylo")) {
    # Call function recursively on all trees in list, then combine results into one data.frame
    return.df <- do.call(rbind, lapply(subtree, sharedMuts))
    return.df$names <- names(subtree)
    return(return.df)
  }

  if (is.null(nu)) {
    nu <- subtree$metadata$nu[1]
    if (is.null(nu)) {
      stop("Need to give a mutation rate (nu) in function call or provide one in params data.frame in tree")
    }
  }

  # Make sure tree is NOT ultrametric
  if (ape::is.ultrametric(subtree)) {
    stop("Tree should be mutation-based, not time-based. Tree should not be ultrametric.")
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

  if (includeStem) {
    message("You have set includeStem = T. Note that we do not include the stem as part of the internal lengths calculation in our work (Johnson et al. 2022)")
  }

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
  descendant_df <- data.frame(
    "Node" = (length(subtree$tip.label) + 2):max(subtree$edge), "Parent" = NA,
    "Edge_length" = NA, "n_cells" = NA
  )


  # Find parent, edge length, and number of descendant cells for each internal node
  for (k in descendant_df$Node) {
    descendants <- subtree$edge[subtree$edge[, 1] == k, 2]
    descendant_df$n_cells[descendant_df$Node == k] <- length(descendants)
    descendant_df$Edge_length[descendant_df$Node == k] <- subtree$edge.length[which(subtree$edge[, 2] == k)]
    descendant_df$Parent[descendant_df$Node == k] <- subtree$edge[which(subtree$edge[, 2] == k), 1]
  }

  # If include Stem is FALSE but tree has a stem, remove stem from calculation
  if (!includeStem & hasStem) {
    descendant_df <- descendant_df[!descendant_df$Parent == stemNode, ]
  }

  # The sum of edge lengths in descendant_df is equal to the total shared mutations
  sharedMutations <- sum(descendant_df$Edge_length)

  # Calculate growth rate and confidence intervals
  growthRate <- n * nu / sharedMutations
  growthRate_lb <- growthRate * (1 + (stats::qnorm(alpha / 2) / sqrt(n)) * (1 + n / sharedMutations))
  growthRate_ub <- growthRate * (1 - (stats::qnorm(alpha / 2) / sqrt(n)) * (1 + n / sharedMutations))

  # Calculate total private (singleton) mutations
  privateMuts <- sum(subtree$edge.length[subtree$edge[, 2] %in% c(1:length(subtree$tip.label))])

  # Check ratio of private to shared mutations
  if (privateMuts / sharedMutations <= 3) {
    warning("Private to shared mutations ratio is less than or equal to 3,
            which means shared mutations method may not be applicable.")
  }

  # Get runtime (including all tests)
  runtime <- proc.time() - ptm

  # return data.frame
  return(data.frame(
    "lowerBound" = growthRate_lb, "estimate" = growthRate,
    "upperBound" = growthRate_ub, "nu" = nu,
    "sharedMutations" = sharedMutations,
    "privateMutations" = privateMuts,
    "extIntRatio" = privateMuts / sharedMutations,
    "n" = n, "alpha" = alpha, "hasStem" = hasStem,
    "includeStem" = includeStem, "runtime_s" = runtime[["elapsed"]],
    "method" = "sharedMuts"
  ))
}



#' Growth rate estimate using Maximum Likelihood
#'
#' @description Uses the approximation that coalescence times H_i are equal to a+b*U_i to
#'     find a and b. b is equal to 1/r, where r is the net growth rate.
#'
#' @param subtree An ape tree subset to include only the clone of interest
#' @param alpha Used for calculation of confidence intervals. 1-alpha confidence
#'     intervals used with default of alpha = 0.05 (95 percent confidence intervals)
#'
#' @return A dataframe including the net growth rate estimate, confidence
#'     intervals, and other important details (runtime, n, etc.)
#' @seealso [cloneRate::internalLengths]
#' @export
#'
#' @examples
#' df <- maxLikelihood(cloneRate::exampleUltraTrees[[1]])
#'
maxLikelihood <- function(subtree, alpha = 0.05) {
  ptm <- proc.time()

  # If we have a list of phylo objects instead of a single phylo objects, call recursively
  if (inherits(subtree, "list") & !inherits(subtree, "phylo")) {
    # Call function recursively on all trees in list, then combine results into one data.frame
    return.df <- do.call(rbind, lapply(subtree, maxLikelihood))
    return.df$names <- names(subtree)
    return(return.df)
  }

  # Basic check on input formatting and alpha value
  inputCheck(subtree, alpha)

  # Get coalescence times
  coal_times <- ape::branching.times(subtree)

  # Log-likelihood function using the approximation for T large
  # params[1]=a, params[2]=r=1/b
  LL <- function(params) {
    a <- params[1]
    r <- params[2]
    U <- (coal_times - a) * r # U_i = (H_i-a)*r
    sigmoid <- 1 / (1 + exp(-U))
    ll <- sum(log(sigmoid)) + sum(log(1 - sigmoid)) + log(r) * length(U)
  }

  # Calculate growth rate by maximizing log likelihood (using maxLik package)
  growthRate <- maxLik::maxLik(LL, start = c(mean(coal_times), .1), )$estimate[2]

  # Get other tree info (lengths)
  extLen <- sum(subtree$edge.length[subtree$edge[, 2] %in% c(1:length(subtree$tip.label))])
  intLen <- suppressWarnings(cloneRate::internalLengths(subtree, includeStem = F)$sumInternalLengths)
  n <- length(subtree$tip.label)
  nodes <- subtree$edge[subtree$edge > n]
  if (1 %in% table(nodes)) {
    hasStem <- T
  } else {
    hasStem <- F
  }

  # Check ratio of external to internal lengths
  if (extLen / intLen <= 3) {
    warning("External to internal lengths ratio is less than or equal to 3,
            which means internal lengths method may not be applicable.")
  }

  # Calculate 1-alpha confidence intervals
  c <- 3 / sqrt(3 + pi^2)
  growthRate_lb <- growthRate * (1 + c * stats::qnorm(alpha / 2) / sqrt(n))
  growthRate_ub <- growthRate * (1 - c * stats::qnorm(alpha / 2) / sqrt(n))

  runtime <- proc.time() - ptm

  return(data.frame(
    "lowerBound" = growthRate_lb, "estimate" = growthRate,
    "upperBound" = growthRate_ub, "sumInternalLengths" = intLen,
    "sumExternalLengths" = extLen, extIntRatio = extLen / intLen,
    "n" = n, "alpha" = alpha, "hasStem" = hasStem,
    "includeStem" = F, "runtime_s" = runtime[["elapsed"]],
    "method" = "maxLike"
  ))
}


#' Check the inputs to growth rate functions
#'
#' @description Check the validity of inputs to growth rate fns, making sure the tree is an
#'     ultrametric ape phylo object with reasonable alpha
#'
#' @param subtree An ape tree subset to include only the clone of interest
#' @param alpha Used for calculation of confidence intervals. 1-alpha confidence
#'     intervals used with default of alpha = 0.05 (95 percent confidence intervals)
#' @keywords internal
#' @return NULL
#'
inputCheck <- function(subtree, alpha) {
  # Must be of class phylo
  if (!inherits(subtree, "phylo")) {
    stop("Tree must be of class phylo. Use as.phylo function to convert if the
    formatting is correct. Otherwise, see ape package documentation
    https://cran.r-project.org/web/packages/ape/ape.pdf")
  }

  # Only works for ultrametric trees
  if (!ape::is.ultrametric(subtree)) {
    stop("Tree is not ultrametric. internalLengths, and maxLike fns.
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






#' Growth rate estimate using Method of Moments (NOT EXPORTED CURRENTLY)
#'
#' @description Provides an estimate for the net growth rate of the clone with confidence
#'     bounds using the method of moments.
#'
#' @param subtree An ape tree subset to include only the clone of interest
#' @param alpha Used for calculation of confidence intervals. 1-alpha confidence
#'     intervals used with default of alpha = 0.05 (95 percent confidence intervals)
#'
#' @returns A dataframe including the net growth rate estimate, confidence
#'     intervals, and other important details (runtime, n, etc.)
#' @seealso [cloneRate::internalLengths()], [cloneRate::maxLikelihood()]
#' @noRd
#' @examples
#' df <- moments(cloneRate::exampleUltraTrees[[1]])
moments <- function(subtree, alpha = 0.05) {
  ptm <- proc.time()

  # If we have a list of phylo objects instead of a single phylo objects, call recursively
  if (inherits(subtree, "list") & !inherits(subtree, "phylo")) {
    # Call function recursively on all trees in list, then combine results into one data.frame
    return.df <- do.call(rbind, lapply(subtree, moments))
    return.df$names <- names(subtree)
    return(return.df)
  }

  # Check if tree has stem
  n <- length(subtree$tip.label)

  # Basic check on input formatting and alpha value
  inputCheck(subtree, alpha)

  # Calculate the growth rate
  growthRate <- (pi / sqrt(3)) * 1 / (stats::sd(ape::branching.times(subtree)))
  growthRate_lb <- growthRate * suppressWarnings(sqrt(1 + 4 * stats::qnorm(alpha / 2) / sqrt(5 * n)))
  if (growthRate_lb == "NaN") {
    growthRate_lb <- 0
  }
  growthRate_ub <- growthRate * sqrt(1 - 4 * stats::qnorm(alpha / 2) / sqrt(5 * n))

  # Get other tree info (lengths)
  extLen <- sum(subtree$edge.length[subtree$edge[, 2] %in% c(1:length(subtree$tip.label))])
  intLen <- suppressWarnings(cloneRate::internalLengths(subtree, includeStem = F)$sumInternalLengths)
  n <- length(subtree$tip.label)
  nodes <- subtree$edge[subtree$edge > n]
  if (1 %in% table(nodes)) {
    hasStem <- T
  } else {
    hasStem <- F
  }

  # Check ratio of external to internal lengths
  if (extLen / intLen <= 3) {
    warning("External to internal lengths ratio is less than or equal to 3,
            which means internal lengths method may not be applicable.")
  }

  runtime <- proc.time() - ptm

  return(data.frame(
    "lowerBound" = growthRate_lb, "estimate" = growthRate,
    "upperBound" = growthRate_ub, "sumInternalLengths" = intLen,
    "sumExternalLengths" = extLen, extIntRatio = extLen / intLen,
    "n" = n, "alpha" = alpha, "hasStem" = hasStem,
    "includeStem" = F, "runtime_s" = runtime[["elapsed"]],
    "method" = "moments"
  ))
}
