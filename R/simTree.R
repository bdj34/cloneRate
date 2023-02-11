#' Simulate ultrametric birth and death branching trees
#'
#' @description Generates a sampled tree (or a list of many sampled trees) from
#'   a supercritial (birth rate > death rate) birth and death branching process
#'   according to the coalescent point process described in "Lambert, A. The
#'   coalescent of a sample from a binary branching process. (2018)."
#'
#' @param a Birth rate
#' @param b Death rate
#' @param cloneAge Clone age. Make sure it's same time units as birth and death
#'   rates
#' @param n Number of samples/tips of the tree to be returned
#' @param precBits Rmpfr param for handling high precision numbers
#' @param addStem Boolean indicating whether to add stem to tree preceding first
#'   split/coalescence
#' @param nTrees Integer indicating the number of trees to generate. Default is
#'   1
#' @param nCores Integer indicating the number of cores to use if parallel pkg
#'  is installed
#'
#' @return An ape object of class "phylo" representing the ultrametric
#'   phylogenetic tree, with elapsed runtime in seconds included. If nTrees
#'   param is greater than 1, a list of objects of class "phylo"
#' @export
#' @importFrom Rmpfr "mpfr"
#' @examples
#' # Generate a single tree
#' tree <- simTree(a = 1, b = 0.5, cloneAge = 20, n = 50)
#'
#' # Generate a list of trees
#' tree_list <- simTree(a = 1, b = 0.5, cloneAge = 20, n = 50, nTrees = 100)
#'
simTree <- function(a, b, cloneAge, n, precBits = 1000, addStem = T,
                    nTrees = 1, nCores = 1) {
  # Store runtime for each tree
  ptm <- proc.time()

  # Check inputs to make sure they make sense
  inputCheck_simTree(
    a = a, b = b, cloneAge = cloneAge, n = n,
    precBits = precBits, addStem = addStem, nTrees = nTrees, nCores = nCores
  )

  # Call recursively to generate nTrees if nTrees > 1
  if (nTrees > 1) {
    # Parallelize if "parallel" pkg avail. (user must set nCores > 1 explicitly)
    if (requireNamespace("parallel", quietly = TRUE)) {
      return.list <- parallel::mclapply(rep(a, nTrees), simTree,
        b = b,
        cloneAge = cloneAge, n = n, precBits = precBits,
        addStem = addStem, nTrees = 1, nCores = 1,
        mc.cores = nCores
      )
    } else {
      return.list <- lapply(rep(a, nTrees), simTree,
        b = b, cloneAge = cloneAge,
        n = n, precBits = precBits, addStem = addStem,
        nTrees = 1
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
  coal_times_mpfr <- vector(length = n - 1)
  for (j in 1:(n - 1)) {
    coal_times_mpfr[j] <- inv_cdf_coal_times(y_mpfr, net_mpfr, a_mpfr, alpha_mpfr, precBits)
  }

  # Convert back to normal numeric (no longer need high precision)
  coal_times <- suppressWarnings(sapply(coal_times_mpfr, Rmpfr::asNumeric))

  # Generate the coalescence intervals for input to ape function
  coal_sorted <- sort(coal_times)
  coal_intervals <- c()
  for (j in 1:length(coal_sorted)) {
    if (j == 1) {
      coal_intervals <- c(coal_intervals, coal_sorted[j])
    } else {
      coal_intervals <- c(coal_intervals, coal_sorted[j] - coal_sorted[j - 1])
    }
  }

  # Make tree with ape function rcoal
  tree <- ape::rcoal(Rmpfr::asNumeric(n_mpfr), rooted = T, br = coal_intervals)

  # Sanity checks (coalescence times must match and be less than cloneAge)
  stopifnot(all(coal_times <= cloneAge))
  if (!all(round(coal_times, 4) %in% round(ape::branching.times(tree), 4))) {
    stop("Unexpected error: coalescence times not matching between input to
         ape::rcoal() and output tree using ape::branching.times()")
  }

  # Add stem starting the tree from zero, rooting the tree appropriately
  if (addStem) {
    tree$edge[tree$edge > n] <- tree$edge[tree$edge > n] + 1
    tree$edge <- rbind(c(n + 1, n + 2), tree$edge)
    tree$edge.length <- c(cloneAge - max(coal_times), tree$edge.length)
    tree$Nnode <- tree$Nnode + 1
  }

  # Add runtime elapsed for making the tree. So fast I know!
  tree$runtime <- proc.time()[["elapsed"]] - ptm[["elapsed"]]

  # Return the tree created from the coalescence times drawn from Lambert distribution
  return(tree)
}


#' Draw a coalescence time from the quantile function
#'
#' @description Draw a single coalescence time from the quantile function by
#'  drawing a random uniform value on 0-1 as input. Uses mpfr numbers to account
#'  for high precision required because y is often very close but not equal to 1.
#'
#' @param y Random variable drawn in simTree(), of class "mpfr"
#' @param net Net growth rate, of class "mpfr"
#' @param a birth rate, of class "mpfr"
#' @param alpha value defined in simTree(), of class "mpfr"
#' @param precBits mpfr param for handling high precision numbers
#' @noRd
#' @return An mpfr number equal to the chosen coalescence time
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
#' @description Make sure that all inputs to the simTree function are
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
#' @return NULL
#'
inputCheck_simTree <- function(a, b, cloneAge, n, precBits, addStem, nTrees,
                               nCores) {
  # Check that we have reasonable inputs
  if (!a > b) {
    stop("simTree function generates trees for supercritical birth-death branching
         processes. Birth rate (a) must be greater than death rate (b).")
  }
  if (!a > 0) {
    stop("simTree function generates trees for supercritical birth-death branching
         processes. Birth rate (a) must be greater than 0.")
  }
  if (b < 0) {
    stop("negative death rate (b) doesn't make sense.")
  }
  if (n > exp((a - b) * cloneAge)) {
    warning("Number of samples (n) is greater than expected population size at
         given net growth rate (a-b) and clone age (cloneAge). This more closely
         resembles a critical rather than the supercritical branching process we
         model. Increase net growth rate or cloneAge, or decrease number of
         samples (n).")
  }
  if (round(n) != n | n < 2) {
    stop("Number of samples must be a positive whole number greater than 1")
  }
  if (round(nTrees) != nTrees | nTrees < 1) {
    stop("Number of trees must be a positive whole number greater than 0")
  }
  if (precBits < 2) {
    stop("precBits for mpfr function is given in bits, and must be at least 2
         bits. For example, double precision corresponds to 53 bits. See Rmpfr
         documentation for details on precBits argument (?Rmpfr::mpfr).")
  }
  if (!inherits(addStem, "logical")) {
    stop("addStem must be logical, indicating whether we want to include a
         'stem' or 'trunk'")
  }
  if (round(nCores) != nCores | nCores < 1) {
    stop("Number of cores must be a positive whole number greater than 0.
         nCores=1 indicates no parallelization.")
  }
  if (nCores > 1 & !requireNamespace("parallel", quietly = T)) {
    warning(paste0("nCores set to ", nCores, " but 'parallel' package is not
                   installed, and simTree() will run without parallelization."))
  }
}
