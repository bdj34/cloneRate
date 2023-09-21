#' Growth rate estimate using the sum of internal lengths
#'
#' @description `internalLengths()` provides an estimate for the net growth rate of the clone with confidence bounds, using the internal lengths method.
#'
#' @param tree An ultrametric tree subset to include only the clone of
#' interest. Alternatively, a list with several such trees.
#' @param alpha Used for calculation of confidence intervals. 1-alpha confidence intervals used with default of alpha = 0.05 (95 percent confidence intervals)
#'
#' @returns A dataframe including the net growth rate estimate, the sum of internal lengths and other important details (clone age estimate, runtime, n, etc.)
#' @seealso [cloneRate::maxLikelihood()], [cloneRate::sharedMuts()] for other
#'  growth rate methods.
#' @export
#' @examples
#' internalLengths(cloneRate::exampleUltraTrees[[1]])
#'
internalLengths <- function(tree, alpha = 0.05) {
  ptm <- proc.time()

  # If we have a list of phylo objects instead of a single phylo objects, call recursively
  if (inherits(tree, "list") & !inherits(tree, "phylo")) {
    # Call function recursively on all trees in list, then combine results into one data.frame
    return.df <- do.call(rbind, lapply(tree, internalLengths, alpha = alpha))
    return.df$cloneName_result <- names(tree)
    return(return.df)
  }

  # Perform basic checks on the input tree
  inputCheck(tree, alpha)

  # Get number of tips
  n <- ape::Ntip(tree)

  # Check if tree has stem
  nodes <- tree$edge[tree$edge > n]
  if (1 %in% table(nodes)) {
    hasStem <- TRUE
    stemNode <- as.numeric(names(which(table(nodes) == 1)))
  } else {
    hasStem <- FALSE
  }

  # Get the number of direct descendants from a node, identifying the nodes with > 2
  countChildren <- table(tree$edge[, 1])

  # Check if tree is binary branching
  if (max(countChildren) > 2) {
    # Throw warning to user
    warningMessage <- paste0(
      "Tree is not binary. Birth-death branching trees should be binary,
       but tree resonstruction from data may lead to  3+ descendants from
       a single parent node. Proceed with caution! Input tree has
      ", max(countChildren), " nodes directly descending from a single
       parent node. A binary tree would only have 2 descendant nodes
       from each parent node."
    )

    if (!is.null(tree$metadata$cloneName_meta)) {
      warningMessage <- paste0(warningMessage, " Tree throwing warning is ", tree$metadata$cloneName_meta[1])
    }
    warning(paste0(warningMessage, "\n"))
  }

  # Get list of descendants from each internal node
  descendant_df <- data.frame(
    "Node" = (n + 2):max(tree$edge), "Parent" = NA,
    "Edge_length" = NA
  )


  # Find parent and edge length preceding each internal node
  for (k in descendant_df$Node) {
    descendant_df$Edge_length[descendant_df$Node == k] <- tree$edge.length[which(tree$edge[, 2] == k)]
    descendant_df$Parent[descendant_df$Node == k] <- tree$edge[which(tree$edge[, 2] == k), 1]
  }

  # If tree has a stem, remove stem from calculation
  if (hasStem) {
    descendant_df <- descendant_df[!descendant_df$Parent == stemNode, ]
  }

  # The sum of edge lengths in descendant_df is equal to the total internal lengths
  intLen <- sum(descendant_df$Edge_length)

  # Calculate growth rate and confidence intervals
  growthRate <- n / intLen
  growthRate_lb <- growthRate * (1 + stats::qnorm(alpha / 2) / sqrt(n))
  growthRate_ub <- growthRate * (1 - stats::qnorm(alpha / 2) / sqrt(n))

  # Calculate total external lengths
  extLen <- sum(tree$edge.length[tree$edge[, 2] %in% c(1:n)])

  # Check ratio of external to internal lengths
  if (extLen / intLen <= 3) {
    warning("External to internal lengths ratio is less than or equal to 3,
            which means internal lengths method may not be applicable. Consider
            using birthDeathMCMC() function, which avoids this issue.\n")
  }

  # Estimate clone age. If tree has stem, take tree age, otherwise estimate by adding 1/r
  if (hasStem) {
    cloneAgeEstimate <- max(ape::branching.times(tree))
  } else {
    cloneAgeEstimate <- max(ape::branching.times(tree)) + 1 / growthRate
  }

  # Get runtime (including all tests)
  runtime <- proc.time() - ptm

  result.df <- data.frame(
    "lowerBound" = growthRate_lb, "estimate" = growthRate,
    "upperBound" = growthRate_ub, "cloneAgeEstimate" = cloneAgeEstimate,
    "sumInternalLengths" = intLen,
    "sumExternalLengths" = extLen, extIntRatio = extLen / intLen,
    "n" = n, "alpha" = alpha, "runtime_s" = runtime[["elapsed"]],
    "method" = "lengths"
  )

  return(result.df)
}



#' Growth rate estimate using the sum of shared mutations assuming a mutation tree
#'
#' @description `sharedMuts()` provides an estimate for the net growth rate of the clone with confidence bounds, using the shared mutations method.
#'
#' @param tree A non-ultrametric ape tree subset to include only the clone of interest
#' @param nu The mutation rate. If none given, sharedMuts() will first look for a `nu` column in a `metadata` data.frame of the tree, and then look for a `nu` in the tree itself. Will throw error if no `nu` given or found.
#' @param alpha Used for calculation of confidence intervals. 1-alpha confidence intervals used with default of alpha = 0.05 (95 percent confidence intervals)
#'
#' @returns A dataframe including the net growth rate estimate, the sum of internal lengths and other important details (clone age estimate, runtime, n, etc.)
#' @seealso [cloneRate::internalLengths()] which is the ultrametric/time-based analogue
#' @export
#' @examples
#' sharedMuts(cloneRate::exampleMutTrees[[1]])
#'
sharedMuts <- function(tree, nu = NULL, alpha = 0.05) {
  ptm <- proc.time()

  # If we have a list of phylo objects instead of a single phylo objects, call recursively
  if (inherits(tree, "list") & !inherits(tree, "phylo")) {
    # Call function recursively on all trees in list, then combine results into one data.frame
    return.df <- do.call(rbind, lapply(tree, sharedMuts))
    return.df$cloneName_result <- names(tree)
    return(return.df)
  }

  # Try to find nu either in metadata data.frame or tree itself
  if (is.null(nu)) {
    nu <- tree$metadata$nu[1]
    if (is.null(nu)) {
      nu <- tree$nu[1]
      if (is.null(nu)) {
        stop("Need to give a mutation rate (nu) in function call or provide one in params data.frame in tree")
      }
    }
  }

  # Make sure tree is NOT ultrametric
  if (ape::is.ultrametric(tree)) {
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

  # Get number of tips
  n <- ape::Ntip(tree)

  # Check if tree has stem
  nodes <- tree$edge[tree$edge > n]
  if (1 %in% table(nodes)) {
    hasStem <- TRUE
    stemNode <- as.numeric(names(which(table(nodes) == 1)))
  } else {
    hasStem <- FALSE
  }

  # Get the number of direct descendants from a node, identifying the nodes with > 2
  countChildren <- table(tree$edge[, 1])

  # Check if tree is binary branching
  if (max(countChildren) > 2) {
    # Throw warning to user
    warningMessage <- paste0(
      "Tree is not binary. Birth-death branching trees should be binary,
       but tree resonstruction from data may lead to  3+ descendants from
       a single parent node. Proceed with caution! Input tree has
      ", max(countChildren), " nodes directly descending from a single
       parent node. A binary tree would only have 2 descendant nodes
       from each parent node."
    )

    if (!is.null(tree$metadata$cloneName_meta)) {
      warningMessage <- paste0(warningMessage, " Tree throwing warning is ", tree$metadata$cloneName_meta[1])
    }
    warning(paste0(warningMessage, "\n"))
  }

  # Get list of descendants from each internal node
  descendant_df <- data.frame(
    "Node" = (n + 2):max(tree$edge), "Parent" = NA,
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
  if (hasStem) {
    descendant_df <- descendant_df[!descendant_df$Parent == stemNode, ]
  }

  # The sum of edge lengths in descendant_df is equal to the total shared mutations
  sharedMutations <- sum(descendant_df$Edge_length)

  # Calculate growth rate and confidence intervals
  growthRate <- n * nu / sharedMutations
  growthRate_lb <- growthRate * (1 + (stats::qnorm(alpha / 2) / sqrt(n)) * (1 + n / sharedMutations))
  growthRate_ub <- growthRate * (1 - (stats::qnorm(alpha / 2) / sqrt(n)) * (1 + n / sharedMutations))

  # Calculate total private (singleton) mutations
  privateMuts <- sum(tree$edge.length[tree$edge[, 2] %in% c(1:n)])

  # Check ratio of private to shared mutations
  if (privateMuts / sharedMutations <= 3) {
    warning("Private to shared mutations ratio is less than or equal to 3,
            which means shared mutations method may not be applicable. If you
            convert the mutation-based tree to a time-based ultrametric tree,
            the birthDeathMCMC() function can be used.\n")
  }

  # Estimate clone age. If tree has stem, take tree age, otherwise estimate by adding 1/r
  if (hasStem) {
    cloneAgeEstimate <- max(ape::branching.times(tree)) / nu
  } else {
    cloneAgeEstimate <- max(ape::branching.times(tree)) / nu + 1 / growthRate
  }

  # Get runtime (including all tests)
  runtime <- proc.time() - ptm

  # return data.frame
  result.df <- data.frame(
    "lowerBound" = growthRate_lb, "estimate" = growthRate,
    "upperBound" = growthRate_ub, "nu" = nu,
    "cloneAgeEstimate" = cloneAgeEstimate,
    "sharedMutations" = sharedMutations,
    "privateMutations" = privateMuts,
    "extIntRatio" = privateMuts / sharedMutations,
    "n" = n, "alpha" = alpha, "runtime_s" = runtime[["elapsed"]],
    "method" = "sharedMuts"
  )

  return(result.df)
}



#' Growth rate estimate using Maximum Likelihood
#'
#' @description Uses the approximation that coalescence times H_i are equal to a+b*U_i to
#'     find a and b. b is equal to 1/r, where r is the net growth rate.
#'
#' @inheritParams internalLengths
#'
#' @returns A dataframe including the net growth rate estimate, confidence
#'     intervals, and other important details (clone age estimate, runtime, n, etc.)
#' @seealso [cloneRate::internalLengths] which uses an alternatvie method for
#'  growth rate estimation from an ultrametric tree.
#' @export
#'
#' @examples
#' df <- maxLikelihood(cloneRate::exampleUltraTrees[[1]])
#'
maxLikelihood <- function(tree, alpha = 0.05) {
  ptm <- proc.time()

  # If we have a list of phylo objects instead of a single phylo object, call recursively
  if (inherits(tree, "list") & !inherits(tree, "phylo")) {
    # Call function recursively on all trees in list, then combine results into one data.frame
    return.df <- do.call(rbind, lapply(tree, maxLikelihood, alpha = alpha))
    return.df$cloneName_result <- names(tree)
    return(return.df)
  }

  # Basic check on input formatting and alpha value
  inputCheck(tree, alpha)

  # Get number of tips
  n <- ape::Ntip(tree)

  # Get the number of direct descendants from a node, identifying the nodes with > 2
  countChildren <- table(tree$edge[, 1])

  # Check if tree is binary branching
  if (!max(countChildren) == 2) {
    # Throw warning to user
    warningMessage <- paste0(
      "Tree is not binary. Birth-death branching trees should be binary,
       but tree resonstruction from data may lead to  3+ descendants from
       a single parent node. Converting tree to binary, but PROCEED WITH CAUTION!
       Input tree has ", max(countChildren), " nodes directly descending
       from a single parent node. A binary tree would only have 2 descendant
       nodes from each parent node."
    )

    if (!is.null(tree$metadata$cloneName_meta)) {
      warningMessage <- paste0(warningMessage, " Tree throwing warning is ", tree$metadata$cloneName_meta)
    }
    warning(paste0(warningMessage, "\n"))

    # Convert tree to binary
    tree <- ape::multi2di(tree)
  }

  # Only take n-1 coal times to avoid using "branching time" from stem node
  coal_times <- sort(ape::branching.times(tree))[c(1:(n - 1))]

  # Negative Log-likelihood function using the approximation for T large
  # params[1]=a, params[2]=r=1/b
  nLL <- function(params) {
    a <- params[1]
    r <- params[2]
    U <- (coal_times - a) * r
    sigmoid <- 1 / (1 + exp(-U))
    ll <- sum(log(sigmoid)) + sum(log(1 - sigmoid)) + log(r) * length(U)
    return(-ll)
  }

  # Calculate growth rate by maximizing log likelihood
  growthRate <- suppressWarnings(stats::optim(
    par = c(mean(coal_times), (pi / sqrt(3)) / stats::sd(coal_times)),
    fn = nLL,
    method = "Nelder-Mead"
  ))$par[2]

  # Get other tree info (lengths)
  extLen <- sum(tree$edge.length[tree$edge[, 2] %in% c(1:n)])
  intLen <- suppressWarnings(cloneRate::internalLengths(tree)$sumInternalLengths)
  nodes <- tree$edge[tree$edge > n]
  if (1 %in% table(nodes)) {
    hasStem <- TRUE
  } else {
    hasStem <- FALSE
  }

  # Check ratio of external to internal lengths
  if (extLen / intLen <= 3) {
    warning("External to internal lengths ratio is less than or equal to 3,
            which means max. likelihood method may not be applicable. Consider
            using birthDeathMCMC() function, which avoids this issue.\n")
  }

  # Calculate 1-alpha confidence intervals
  c <- 3 / sqrt(3 + pi^2)
  growthRate_lb <- growthRate * (1 + c * stats::qnorm(alpha / 2) / sqrt(n))
  growthRate_ub <- growthRate * (1 - c * stats::qnorm(alpha / 2) / sqrt(n))

  # Estimate clone age. If tree has stem, take tree age, otherwise estimate by adding 1/r
  if (hasStem) {
    cloneAgeEstimate <- max(ape::branching.times(tree))
  } else {
    cloneAgeEstimate <- max(ape::branching.times(tree)) + 1 / growthRate
  }

  runtime <- proc.time() - ptm

  return(data.frame(
    "lowerBound" = growthRate_lb, "estimate" = growthRate,
    "upperBound" = growthRate_ub, "cloneAgeEstimate" = cloneAgeEstimate,
    "sumInternalLengths" = intLen,
    "sumExternalLengths" = extLen, extIntRatio = extLen / intLen,
    "n" = n, "alpha" = alpha, "runtime_s" = runtime[["elapsed"]],
    "method" = "maxLike"
  ))
}





#' Growth rate estimate using MCMC
#'
#' @description Uses Rstan and the No U-turn sampler to approximate the
#'  growth rate using the likelihood from Stadler 2009 "On incomplete sampling
#'  under birth–death models and connections to the sampling-based coalescent"
#'
#' @inheritParams internalLengths
#' @param maxGrowthRate Sets upper bound on birth rate. Default is 4 but this
#'  will depend on the nature of the data
#' @param verbose TRUE or FALSE, should the Rstan MCMC intermediate output and progress be printed?
#' @param nChains Number of chains to run in MCMC. Default is 4
#' @param nCores Number of cores to perform MCMC. Default is 1, but chains can
#'  be run in parallel
#' @param chainLength Number of iterations for each chain in MCMC. Default is
#'  2000 (1000 warm-up + 1000 sampling), increase if stan tells you to
#'
#' @returns A dataframe including the net growth rate estimate, confidence
#'     intervals, and other important details (clone age estimate, runtime, n,
#'     etc.)
#' @seealso [cloneRate::internalLengths] [cloneRate::maxLikelihood] which use
#'  alternative methods for growth rate estimation from an ultrametric tree.
#' @export
#'
#' @examples
#' \donttest{
#' df <- birthDeathMCMC(cloneRate::exampleUltraTrees[[1]])
#' }
#'
birthDeathMCMC <- function(tree, maxGrowthRate = 4, alpha = 0.05,
                           verbose = TRUE, nChains = 4,
                           nCores = 1, chainLength = 2000) {
  # If we have a list of phylo objects instead of a single phylo object, call recursively
  if (inherits(tree, "list") & inherits(tree[[1]], "phylo")) {
    # Run birth-death MCMC model many times. Parallelize if possible
    if (requireNamespace("parallel")) {
      df <- do.call(rbind, parallel::mclapply(tree, runStan,
        stanModel = stanmodels$bdSampler,
        maxGrowthRate = maxGrowthRate, alpha = alpha,
        verbose = verbose, nChains = nChains,
        nCores = if (nChains == nCores) {
          nCores
        } else {
          1
        },
        chainLength = chainLength,
        mc.cores = if (nChains == nCores) {
          1
        } else {
          nCores
        }
      ))
      df$cloneName_result <- names(tree)
      return(df)
    } else {
      df <- do.call(rbind, lapply(tree, runStan,
        stanModel = stanmodels$bdSampler,
        maxGrowthRate = maxGrowthRate, alpha = alpha,
        verbose = verbose, nChains = nChains,
        nCores = nCores, chainLength = chainLength
      ))
      df$cloneName_result <- names(tree)
      return(df)
    }
  } else {
    # Run stan model once
    df <- runStan(tree,
      stanModel = stanmodels$bdSampler,
      maxGrowthRate = maxGrowthRate, alpha = alpha,
      verbose = verbose, nChains = nChains,
      nCores = nCores, chainLength = chainLength
    )
    return(df)
  }
}


#' Run stan model
#' @noRd
#'
#' @description Takes a compiled stan model and runs it. Same params as
#'  birthDeathMCMC and one additional param, stanModel, which is the
#'  object of class stanmodel (essentially compiled stan). See rstan package
#'  for more details
#'
#' @inheritParams birthDeathMCMC
#' @param stanModel Compiled stan model, generated using rstan::stan_model
#'
#' @keywords internal
#' @returns A dataframe including the net growth rate estimate, confidence
#'  intervals, and other important details (clone age estimate, runtime, n,
#'  etc.)
#'
#'
runStan <- function(tree, stanModel, maxGrowthRate = 4, alpha = 0.05,
                    verbose = TRUE, nChains = 3,
                    nCores = 1, chainLength = 2000) {
  # Keep track of time to run each instance (excluding compile time)
  ptm <- proc.time()

  # Basic check on input formatting and alpha value
  inputCheck(tree, alpha)

  # Get number of tips
  n <- ape::Ntip(tree)

  # Get the number of direct descendants from a node, identifying the nodes with > 2
  countChildren <- table(tree$edge[, 1])

  # Check if tree is binary branching
  if (!max(countChildren) == 2) {
    # Throw warning to user
    warningMessage <- paste0(
      "Tree is not binary. Birth-death branching trees should be binary,
       but tree resonstruction from data may lead to  3+ descendants from
       a single parent node. Converting tree to binary, but PROCEED WITH CAUTION!
       Input tree has ", max(countChildren), " nodes directly descending
       from a single parent node. A binary tree would only have 2 descendant
       nodes from each parent node."
    )

    if (!is.null(tree$metadata$cloneName_meta)) {
      warningMessage <- paste0(warningMessage, " Tree throwing warning is ", tree$metadata$cloneName_meta)
    }
    warning(paste0(warningMessage, "\n"))

    # Convert tree to binary
    tree <- ape::multi2di(tree)
  }

  # Get internal lengths from our lengths function
  resultLengths <- suppressWarnings(internalLengths(tree))

  # Get n-1 coalescence times, removing the nth if given a tree with a stem
  coal_times <- sort(ape::branching.times(tree))[c(1:(n - 1))]
  nCoal <- length(coal_times)

  # Set input data
  inData <- list(
    "nCoal" = nCoal,
    "t" = coal_times,
    "upperLambda" = maxGrowthRate
  )

  # Run Rstan, setting variable
  if (verbose) {
    stanr <- tryCatch.W.E({
      rstan::sampling(stanModel,
        data = inData,
        chains = nChains,
        cores = nCores,
        iter = chainLength,
        verbose = TRUE
      )
    })
  } else {
    stanr <- tryCatch.W.E({
      rstan::sampling(stanModel,
        data = inData,
        chains = nChains,
        cores = nCores,
        iter = chainLength,
        refresh = 0
      )
    })
  }

  outList <- list(
    posterior = rstan::extract(stanr$value),
    res = stanr$value,
    tree = tree,
    dat = inData
  )

  warningMessage <- paste(unlist(stanr$warning), collapse = " ____ ")

  # Get growth rate and 95% CI, alos rough estimate of sampling probability
  ptile <- c(alpha / 2, 0.5, 1 - alpha / 2)
  growthRateVec <- stats::quantile(outList$posterior$lambda - outList$posterior$mu, ptile)
  rhoVec <- exp(stats::quantile(outList$posterior$lgRho, ptile))

  # Rough estimate of clone age
  cloneAgeEstimate <- max(coal_times) + 1 / growthRateVec[2]

  runtime <- proc.time() - ptm

  return(data.frame(
    "lowerBound" = growthRateVec[1], "estimate" = growthRateVec[2],
    "upperBound" = growthRateVec[3], "cloneAgeEstimate" = cloneAgeEstimate,
    "sumInternalLengths" = resultLengths$sumInternalLengths,
    "sumExternalLengths" = resultLengths$sumExternalLengths,
    "extIntRatio" = resultLengths$extIntRatio,
    "n" = n, "alpha" = alpha, "runtime_s" = runtime[["elapsed"]],
    "method" = "BirthDeathMCMC", "samplingProbBallpark" = rhoVec[2],
    "chainLength" = chainLength, "nChains" = nChains, "nCores" = nCores,
    "warningMessage" = warningMessage
  ))
}




#' Check the inputs to growth rate functions
#'
#' @description Check the validity of inputs to growth rate fns, making sure the tree is an
#'     ultrametric ape phylo object with reasonable alpha
#'
#' @param tree An ape tree subset to include only the clone of interest
#' @param alpha Used for calculation of confidence intervals. 1-alpha confidence
#'     intervals used with default of alpha = 0.05 (95 percent confidence intervals)
#' @keywords internal
#' @returns NULL
#'
inputCheck <- function(tree, alpha) {
  # Must be of class phylo
  if (!inherits(tree, "phylo")) {
    stop("Tree must be of class phylo. Use as.phylo function to convert if the
    formatting is correct. Otherwise, see ape package documentation
    https://cran.r-project.org/web/packages/ape/ape.pdf")
  }

  # Only works for ultrametric trees
  if (!ape::is.ultrametric(tree)) {
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
#' @inheritParams internalLengths
#'
#' @returns A dataframe including the net growth rate estimate, confidence
#'     intervals, and other important details (clone age estimate, runtime, n, etc.)
#' @seealso [cloneRate::internalLengths()], [cloneRate::maxLikelihood()]
#' @noRd
#' @examples
#' df <- moments(cloneRate::exampleUltraTrees[[1]])
moments <- function(tree, alpha = 0.05) {
  ptm <- proc.time()

  # If we have a list of phylo objects instead of a single phylo objects, call recursively
  if (inherits(tree, "list") & !inherits(tree, "phylo")) {
    # Call function recursively on all trees in list, then combine results into one data.frame
    return.df <- do.call(rbind, lapply(tree, moments, alpha = alpha))
    return.df$cloneName_result <- names(tree)
    return(return.df)
  }

  # Get number of tips
  n <- ape::Ntip(tree)

  # Basic check on input formatting and alpha value
  inputCheck(tree, alpha)

  # Get the number of direct descendants from a node, identifying the nodes with > 2
  countChildren <- table(tree$edge[, 1])

  # Check if tree is binary branching
  if (!max(countChildren) == 2) {
    # Throw warning to user
    warningMessage <- paste0(
      "Tree is not binary. Birth-death branching trees should be binary,
       but tree resonstruction from data may lead to  3+ descendants from
       a single parent node. Converting tree to binary, but PROCEED WITH CAUTION!
       Input tree has ", max(countChildren), " nodes directly descending
       from a single parent node. A binary tree would only have 2 descendant
       nodes from each parent node."
    )

    if (!is.null(tree$metadata$cloneName_meta)) {
      warningMessage <- paste0(warningMessage, " Tree throwing warning is ", tree$metadata$cloneName_meta)
    }
    warning(paste0(warningMessage, "\n"))

    # Convert tree to binary
    tree <- ape::multi2di(tree)
  }

  # Only take n-1 coal times to avoid using "branching time" from stem node
  coal_times <- sort(ape::branching.times(tree))[c(1:(n - 1))]


  # Calculate the growth rate
  growthRate <- (pi / sqrt(3)) * 1 / (stats::sd(coal_times))
  growthRate_lb <- growthRate * suppressWarnings(sqrt(1 + 4 * stats::qnorm(alpha / 2) / sqrt(5 * n)))
  if (growthRate_lb == "NaN") {
    growthRate_lb <- 0
  }
  growthRate_ub <- growthRate * sqrt(1 - 4 * stats::qnorm(alpha / 2) / sqrt(5 * n))

  # Get other tree info (lengths)
  extLen <- sum(tree$edge.length[tree$edge[, 2] %in% c(1:n)])
  intLen <- suppressWarnings(cloneRate::internalLengths(tree)$sumInternalLengths)
  nodes <- tree$edge[tree$edge > n]
  if (1 %in% table(nodes)) {
    hasStem <- TRUE
  } else {
    hasStem <- FALSE
  }

  # Check ratio of external to internal lengths
  if (extLen / intLen <= 3) {
    warning("External to internal lengths ratio is less than or equal to 3,
            which means moments method may not be applicable.\n")
  }

  # Estimate clone age. If tree has stem, take tree age, otherwise estimate by adding 1/r
  if (hasStem) {
    cloneAgeEstimate <- max(ape::branching.times(tree))
  } else {
    cloneAgeEstimate <- max(ape::branching.times(tree)) + 1 / growthRate
  }

  runtime <- proc.time() - ptm

  return(data.frame(
    "lowerBound" = growthRate_lb, "estimate" = growthRate,
    "upperBound" = growthRate_ub, "cloneAgeEstimate" = cloneAgeEstimate,
    "sumInternalLengths" = intLen,
    "sumExternalLengths" = extLen, extIntRatio = extLen / intLen,
    "n" = n, "alpha" = alpha, "runtime_s" = runtime[["elapsed"]],
    "method" = "moments"
  ))
}


#' Capture error and warning messages
#'
#' @description R utility function. Run
#'  file.show(system.file("demo/error.catching.R")) for details)
#'
#' @noRd
#' @param expr Expression of R code to run can capture output + warnings from
#' @keywords internal
#'
#' @returns list with value as expression output and warning as warning
#'
tryCatch.W.E <- function(expr) {
  W <- NULL
  w.handler <- function(w) { # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(
    value = withCallingHandlers(tryCatch(expr, error = function(e) e),
      warning = w.handler
    ),
    warning = W
  )
}
