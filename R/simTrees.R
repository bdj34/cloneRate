#' Simulate ultrametric birth and death branching trees
#'
#' @description Generates a sampled tree (or a list of many sampled trees) from
#'   a supercritial (birth rate > death rate) birth and death branching process
#'   according to the coalescent point process described in "Lambert, A. The
#'   coalescent of a sample from a binary branching process. (2018)."
#'
#' @param a Birth rate or vector of birth rates of length 'nTrees'
#' @param b Death rate or vector of death rates of length 'nTrees'
#' @param cloneAge Clone age or vector of clone ages of length 'nTrees'. Make
#'   sure it's same time units as birth and death rates
#' @param n Number of samples/tips of the tree to be returned. Can be a vector
#'   of length 'nTrees' as well.
#' @param nTrees Integer indicating the number of trees to generate. Default is
#'   1
#' @param precBits Rmpfr param for handling high precision numbers. Needed for
#'   drawing the coalescence times. Can be a vector of length 'nTrees', though
#'   it is not recommended
#' @param addStem Boolean indicating whether to add stem to tree preceding first
#'   split/coalescence. Can also be a vector of length 'nTrees'
#' @param nCores Integer indicating the number of cores to use if parallel pkg
#'  is installed. Default is 1.
#'
#' @returns An ape object of class "phylo" representing the ultrametric
#'   phylogenetic tree with edge lengths in units of time. Tree metadata is
#'   located in the 'metadata' data.frame included in each "phylo" object. If
#'   'nTrees' param is greater than 1, simUltra returns a list of objects of
#'   such objects of class "phylo".
#' @export
#' @importFrom Rmpfr "mpfr"
#' @examples
#' # Generate a single tree
#' tree <- simUltra(a = 1, b = 0.5, cloneAge = 20, n = 50)
#'
#' # Generate a list of trees
#' tree_list <- simUltra(a = 1, b = 0.5, cloneAge = 20, n = 50, nTrees = 3)
#'
simUltra <- function(a, b, cloneAge, n, nTrees = 1,
                     precBits = 1000, addStem = FALSE, nCores = 1) {
  # Store runtime for each tree
  ptm <- proc.time()

  # Set nTrees equal to "SKIP_TESTS" on recursive runs to skip repeated checks
  if (nTrees == "SKIP_TESTS") {
    nTrees <- 1
  } else {
    # Make sure length of params is either one or equal to 'nTrees'
    if (!all(unlist(lapply(list(a, b, cloneAge, n), length)) == 1 |
      unlist(lapply(list(a, b, cloneAge, n), length)) == nTrees)) {
      stop(paste0("Input parameters must be length 1 or length equal to the value
                of param 'nTrees', which is ", nTrees))
    }

    # Check inputs to make sure they make sense
    inputCheck_simTree(
      a = a, b = b, cloneAge = cloneAge, n = n,
      precBits = precBits, addStem = addStem, nTrees = nTrees, nCores = nCores
    )

    # If a is length 1, make it a vector of length nTrees for first argument to mapply
    if (length(a) == 1) {
      a <- rep(a, nTrees)
    } else {
      a <- a
    }
  }

  # Call recursively to generate nTrees if nTrees > 1
  if (nTrees > 1) {
    # Parallelize if "parallel" pkg avail. (user must set nCores > 1 explicitly)
    if (requireNamespace("parallel", quietly = TRUE)) {
      return.list <- parallel::mcmapply(simUltra,
        a = a,
        b = b,
        cloneAge = cloneAge, n = n, precBits = precBits,
        addStem = addStem, nTrees = "SKIP_TESTS", nCores = 1,
        mc.cores = nCores, SIMPLIFY = FALSE
      )
    } else {
      return.list <- mapply(
        FUN = simUltra, a = a,
        b = b, cloneAge = cloneAge,
        n = n, precBits = precBits, addStem = addStem,
        nTrees = "SKIP_TESTS", SIMPLIFY = FALSE
      )
    }
    return(return.list)
  }

  # Convert params to high precision mpfr
  cloneAge_mpfr <- mpfr(cloneAge, precBits)
  a_mpfr <- mpfr(a, precBits)
  b_mpfr <- mpfr(b, precBits)
  net_mpfr <- a_mpfr - b_mpfr
  n_mpfr <- mpfr(n, precBits)
  one <- mpfr(1, precBits)

  # Define alpha and ensure it doesn't map to 1 at given precision
  alpha_mpfr <- (a_mpfr * (exp(net_mpfr * cloneAge_mpfr) - one)) / (a_mpfr * exp(net_mpfr * cloneAge_mpfr) - b_mpfr)
  if (alpha_mpfr == one) {
    stop("alpha value is equal to 1 due to insufficient machine precision. This
          will lead to NaN coalescence times. Increase param 'precBits'.")
  }

  # Draw Y = y from the inverse CDF
  uniform_rv <- mpfr(stats::runif(1, min = 0, max = 1), precBits)
  y_mpfr <- ((one - alpha_mpfr) * (uniform_rv**(one / n_mpfr))) / (one - alpha_mpfr * uniform_rv**(one / n_mpfr))

  # Generate the coalescence times
  coal_times_mpfr <- sapply(rep(y_mpfr, n - 1), inv_cdf_coal_times,
    net = net_mpfr,
    a = a_mpfr, alpha = alpha_mpfr, precBits = precBits
  )

  # Convert back to normal numeric (no longer need high precision)
  coal_times <- suppressWarnings(sapply(coal_times_mpfr, Rmpfr::asNumeric))

  # Convert coal times into tree by randomly merging lineages
  tree <- coal_to_tree(coal_times)

  # Add stem starting the tree from zero, rooting the tree appropriately
  if (addStem) {
    tree$edge[tree$edge > n] <- tree$edge[tree$edge > n] + 1
    tree$edge <- rbind(c(n + 1, n + 2), tree$edge)
    tree$edge.length <- c(cloneAge - max(coal_times), tree$edge.length)
    tree$Nnode <- tree$Nnode + 1
  }

  # Sanity checks (coalescence times must match and be less than cloneAge)
  stopifnot(all(coal_times <= cloneAge))
  if (!all(round(coal_times, 4) %in% round(ape::branching.times(tree), 4))) {
    stop("Unexpected error: coalescence times not matching between drawn times
         and output tree using ape::branching.times()")
  }

  # Add metadata for making the tree
  runtime <- proc.time()[["elapsed"]] - ptm[["elapsed"]]
  tree$metadata <- data.frame(
    "r" = a - b, "a" = a, "b" = b, "cloneAge" = cloneAge,
    "n" = n, "runtime_seconds" = runtime, "addStem" = addStem
  )

  # Return the tree created from the coalescence times drawn from Lambert distribution
  return(tree)
}






#' Simulate mutation-based birth and death branching trees
#'
#' @description Generates a sampled tree (or a list of many sampled trees) from
#'   a supercritial (birth rate > death rate) birth and death branching process
#'   according to the coalescent point process described in "Lambert, A. The
#'   coalescent of a sample from a binary branching process. (2018)." Edge
#'   lengths will be in units of mutations, assuming poissonian mutation
#'   accumulation. Essentially a wrapper combining simUltra() and ultra2mut()
#'   functions into one step.
#'
#' @param a Birth rate
#' @param b Death rate
#' @param cloneAge Clone age. Make sure it's same time units as birth and death
#'   rates
#' @param n Number of samples/tips of the tree to be returned
#' @param nu Mutation rate in units of mutations per unit time. Make sure time
#'   units are consistent with birth and death rates and cloneAge
#' @param nTrees Integer indicating the number of trees to generate. Default is
#'   1.
#' @param precBits Rmpfr param for handling high precision numbers. Needed to
#'   draw coalescence times.
#' @param addStem Boolean indicating whether to add stem to tree preceding first
#'   split/coalescence
#' @param nCores Integer indicating the number of cores to use if parallel pkg
#'  is installed. Default is 1.
#'
#' @returns An ape object of class "phylo" representing the ultrametric
#'   phylogenetic tree with edge lengths in units of time. Tree metadata is
#'   located in the 'metadata' data.frame included in each "phylo" object. If
#'   'nTrees' param is greater than 1, simUltra returns a list of objects of
#'   such objects of class "phylo".
#' @export
#' @importFrom Rmpfr "mpfr"
#' @examples
#' # Generate a single mutation-based tree with a specified mutation rate
#' tree <- simMut(a = 1, b = 0.5, cloneAge = 40, n = 50, nu = 10)
#'
#' # Generate a list of mutation-based trees with a range of mutation rates
#' tree_list <- simMut(
#'   a = 1, b = 0.5, cloneAge = 40, n = 50,
#'   nu = stats::runif(n = 3, min = 10, max = 20), nTrees = 3
#' )
#'
simMut <- function(a, b, cloneAge, n, nu, nTrees = 1,
                   precBits = 1000, addStem = FALSE, nCores = 1) {
  # Generate ultrametric, time-based trees
  ultraTrees <- simUltra(
    a = a, b = b, cloneAge = cloneAge, n = n,
    precBits = precBits, addStem = addStem, nTrees = nTrees,
    nCores = nCores
  )

  # Convert from ultrametric, time-based trees to mutation based trees
  ultra2mut(ultraTrees, nu = nu)
}




#' Add poissonian mutations to an ultrametric tree(s)
#'
#' @description Takes an ultrametric tree of class "phylo" (or a list of such
#' trees) and draws new edge lengths in units of mutations, with the mean of
#' each new edge length equal to the old edge length multiplied by the mutation
#' rate. Mutation rate can be set or drawn from a uniform distribution.
#'
#' @param tree A single tree or list of trees of class "phylo", with edge
#'   lengths in units of time
#' @param nu Mutation rate in units of mutations per unit time. Can also be a
#'   vector of mutation rates with length equal to the number of input trees.
#'   Make sure time units are consistent in nu and tree$edge.length
#'
#' @returns An ape object of class "phylo" representing the phylogenetic tree
#'   with edge lengths in units of mutations. Value of mutation rate will be
#'   added to 'metadata' data.frame of output tree if such a data.frame exists
#'   in the input tree. Otherwise, mutation rate value will be added to "phylo"
#'   object directly. If input is a list of trees, ultra2mut() will return a
#'   list of such "phylo" objects.
#'
#' @export
#' @examples
#' # Convert the time-based, ultrametric example trees into mutation-based trees
#' mutTrees <- ultra2mut(exampleUltraTrees,
#'   nu = stats::runif(n = length(exampleUltraTrees), min = 10, max = 20)
#' )
#'
ultra2mut <- function(tree, nu) {
  # If we have a list of phylo objects instead of a single phylo objects, call recursively
  if (inherits(tree, "list") & !inherits(tree, "phylo")) {
    if (!(length(nu) == 1 | length(nu) == length(tree))) {
      stop("Length of nu must be 1 or equal to the number of input trees")
    }
    # Call function recursively on all trees in list, then combine results into one data.frame
    return.df <- mapply(ultra2mut, tree, nu = nu, SIMPLIFY = FALSE)
    return(return.df)
  }

  # Make sure tree has edge lengths
  if (is.null(tree$edge.length)) {
    stop("Tree doesn't have edge lengths so there is no way to generate random
         number of mutations. Use a 'phylo' object with an 'edge.length' vector.")
  }

  # Adjust edge lengths of the tree using poisson random draws
  tree$edge.length <- sapply(tree$edge.length, function(x) {
    stats::rpois(1, nu * x)
  })

  # If tree has a 'metadata' data.frame, add 'nu' to it, otherwise add 'nu' to tree directly
  if (is.null(tree$metadata)) {
    tree$nu <- nu
  } else {
    tree$metadata$nu <- nu
  }

  return(tree)
}



#' Generate tree from coalescence times
#'
#' @description generates a tree from a vector of coalescence times by randomly
#'  merging lineages.
#'
#' @param coal_times A numeric vector of coalescence times
#'
#' @returns An ape object of class "phylo" representing the ultrametric
#'   phylogenetic tree with edge lengths in units of time.
#' @export
#'
#' @examples
#' # Generate an ape phylo tree with n tips from a vector of n-1 coalescence times
#' randomCoalTimes <- c(9.3, 7.8, 10.15, 11.23, 9.4, 8.8, 10.01, 13)
#' tree <- coal_to_tree(randomCoalTimes)
#'
coal_to_tree <- function(coal_times) {
  # coal_times must be a vector of numbers
  if (!inherits(coal_times, "numeric") | length(coal_times) < 2) {
    stop("coal_times input to coal_to_tree() function must a numeric vector")
  }

  # Get number of tips, n
  n <- length(coal_times) + 1

  # Set number of total edges and initialize edge and edge.length
  numEdges <- 2 * n - 2
  edge.length <- rep(0, numEdges)
  edge <- matrix(NA, nrow = numEdges, ncol = 2)

  # Fill heights vec to keep track of height of nodes
  heights <- rep(0, numEdges - 1)

  # Sort coalescence times (smallest to largest)
  coal_times_sorted <- sort(coal_times, decreasing = FALSE)
  possibleChildren <- as.integer(c(1:n))
  currentNode <- as.integer(2 * n - 1)

  # Loop through n-1 internal nodes
  for (i in 1:(n - 1)) {
    # Sample the children
    children <- sample(possibleChildren, size = 2, replace = FALSE)

    # Go to next open row
    row <- which(is.na(edge[, 1]))[c(1, 2)]

    # Fill second column with children and first with node
    edge[row, 2] <- children
    edge[row, 1] <- currentNode

    # Set edge.length as diff. between coal time and children height
    edge.length[row] <- coal_times_sorted[i] - heights[children]

    # Set height of current node
    heights[currentNode] <- coal_times_sorted[i]

    # Add current node to list of possible children
    possibleChildren <- c(possibleChildren[!possibleChildren %in% children], currentNode)

    # Move on to the next current node
    currentNode <- currentNode - 1L
  }

  # Make tree as list
  tree <- list(edge = edge, edge.length = edge.length, Nnode = as.integer(n - 1))
  tree$tip.label <- sample(paste0("t", c(1:n)), replace = FALSE)

  # Set class
  class(tree) <- "phylo"

  return(tree)
}



#' Draw a coalescence time from the quantile function
#'
#' @description Draw a single coalescence time from the quantile function by
#'  drawing a random uniform value on 0-1 as input. Uses mpfr numbers to account
#'  for high precision required because y is often very close but not equal to 1.
#'
#' @param y Random variable drawn in simUltra(), of class "mpfr"
#' @param net Net growth rate, of class "mpfr"
#' @param a birth rate, of class "mpfr"
#' @param alpha value defined in simUltra(), of class "mpfr"
#' @param precBits mpfr param for handling high precision numbers
#' @noRd
#' @returns An mpfr number equal to the chosen coalescence time
#' @importFrom Rmpfr "mpfr"
#'
inv_cdf_coal_times <- function(y, net, a, alpha, precBits) {
  one <- mpfr(1, precBits)
  rv <- mpfr(stats::runif(1, min = 0, max = 1), precBits)
  phi <- (alpha * rv) / (a * (one - alpha * (one - y)))
  return((-one / net) * log((one - y * a * phi) / (one + (net - y * a) * phi)))
}





#' Check that inputs are reasonable
#'
#' @description Make sure that all inputs to the simUltra function are
#'  reasonable, making it less likely that the user makes a mistake. Provide
#'  helpful error messages if something is wrong.
#'
#' @param a Birth rate
#' @param b Death rate
#' @param cloneAge Clone age (make sure it's same time units as birth and death
#'   rates)
#' @param n Number of samples (number of tips of the tree to be returned)
#' @param precBits Rmpfr param for handling high precision numbers
#' @param addStem Boolean indicating whether to add stem to tree preceding first
#'   split/coalescence
#' @param nTrees Integer indicating the number of trees to be generated
#' @noRd
#' @returns NULL
#'
inputCheck_simTree <- function(a, b, cloneAge, n, precBits, addStem, nTrees,
                               nCores) {
  # Check that we have reasonable inputs
  if (!all(a > b)) {
    stop("simUltra function generates trees for supercritical birth-death branching
         processes. Birth rate (a) must be greater than death rate (b).")
  }
  if (!all(a > 0)) {
    stop("simUltra function generates trees for supercritical birth-death branching
         processes. Birth rate (a) must be greater than 0.")
  }
  if (any(b < 0)) {
    stop("negative death rate (b) doesn't make sense.")
  }
  if (any(n > exp((a - b) * cloneAge))) {
    warning("Number of samples (n) is greater than expected population size at
         given net growth rate (a-b) and clone age (cloneAge). This more closely
         resembles a critical rather than the supercritical branching process we
         model. Increase net growth rate or cloneAge, or decrease number of
         samples (n).")
  }
  if (any(round(n) != n) | any(n < 2)) {
    stop("Number of samples must be a positive whole number greater than 1")
  }
  if (any(round(nTrees) != nTrees) | any(nTrees < 1) | length(nTrees) != 1) {
    stop("Number of trees must be a positive whole number greater than 0")
  }
  if (any(precBits < 2)) {
    stop("precBits for mpfr function is given in bits, and must be at least 2
         bits. For example, double precision corresponds to 53 bits. See Rmpfr
         documentation for details on precBits argument (?Rmpfr::mpfr).")
  }
  if (!all(inherits(addStem, "logical"))) {
    stop("addStem must be logical, indicating whether we want to include a
         'stem' or 'trunk'")
  }
  if (any(round(nCores) != nCores) | any(nCores < 1)) {
    stop("Number of cores must be a positive whole number greater than 0.
         nCores=1 indicates no parallelization.")
  }
  if (any(nCores > 1) & !all(requireNamespace("parallel", quietly = TRUE))) {
    warning(paste0("nCores set to ", nCores, " but 'parallel' package is not
                   installed, and simUltra() will run without parallelization."))
  }
  if (any(n < 4)) {
    stop("Number of samples, n, must be greater than 3.")
  }
}
