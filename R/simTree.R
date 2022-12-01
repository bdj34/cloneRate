#' @title simTree
#'
#' @description Generates a sampled tree from a supercritial (birth rate >
#'     death rate) birth and death branching process according to the coalescent
#'     point process described in "Lambert, A. The coalescent of a sample from a
#'     binary branching process. (2018)."
#'
#' @param a Birth rate
#' @param b Death rate
#' @param cloneAge Clone age
#' @param n Number of samples (number of tips of the tree to be returned)
#' @param precision Rmpfr param for handling high precision numbers
#'
#' @return An ape tree object
#' @examples
#' tree <- simTree(a=1, b=0.5, cloneAge = 20, n = 50)
#' @export
#' @importFrom Rmpfr "mpfr"
#' @importFrom Rmpfr "asNumeric"
#' @importFrom ape "rcoal"
#' @importFrom ape "branching.times"
simTree <- function(a, b, cloneAge, n, precision = 10000){

  # Check that we have reasonable inputs
  if (! a > b){
    stop("This function generates trees for supercritical birth-death branching
         processes. Birth rate (a) must be greater than death rate (b).")
  }
  # Convert params to high precision mpfr
  cloneAge_mpfr <- mpfr(cloneAge, precision)
  a_mpfr <- mpfr(a, precision)
  b_mpfr <- mpfr(b, precision)
  net_mpfr <- a_mpfr - b_mpfr
  n_mpfr <- mpfr(n, precision)
  one <- mpfr(1, precision)



  #Define alpha
  alpha_mpfr <- (a_mpfr*(exp(net_mpfr*cloneAge_mpfr)-one))/(a_mpfr*exp(net_mpfr*cloneAge_mpfr)-b_mpfr)

  #Define inverse CDF for the coalescence times
  inv_cdf_coal_times <- function(y, net. = net_mpfr, a. = a_mpfr, alpha. = alpha_mpfr){
    rv <- mpfr(runif(1), precision)
    phi <- (alpha.*rv)/(a.*(one-alpha.*(one-y)))
    return((-one/net.)*log((one-y*a.*phi)/(one+(net.-y*a.)*phi)))
  }

  #Draw Y = y from the inverse CDF
  uniform_rv <- mpfr(runif(1), precision)
  y_mpfr <- ((one-alpha_mpfr)*(uniform_rv**(one/n_mpfr)))/(one-alpha_mpfr*uniform_rv**(one/n_mpfr))

  # Generate the coalesence times
  coal_times_mpfr <- vector(length = n-1)
  for (j in 1:(n-1)){
    coal_times_mpfr[j] <- inv_cdf_coal_times(y_mpfr)
  }

  # Convert back to normal numeric (no longer need high precision)
  coal_times <- suppressWarnings(sapply(coal_times_mpfr, asNumeric))

  # Generate the coalescence intervals for input to ape function
  coal_sorted <- sort(coal_times)
  coal_intervals <- c()
  for (j in 1:length(coal_sorted)) {
    if (j == 1){
      coal_intervals <- c(coal_intervals, coal_sorted[j])
    } else {
      coal_intervals <- c(coal_intervals, coal_sorted[j] - coal_sorted[j-1])
    }
  }

  # Make tree with ape function
  tree <- rcoal(asNumeric(n_mpfr), rooted = T, br = coal_intervals)

  # Sanity checks
  stopifnot(all(round(coal_times,4) %in% round(branching.times(tree),4)))
  stopifnot(all(coal_times <= cloneAge))
  #if (! all(round(coal_times,4) %in% round(branching.times(tree),4))) {
  #  print("Unexpected error: coal. times not matching!")
  #  return(NULL)
  #}

  # Add outgroup starting the tree from zero then remove, rooting the tree
  #tree <- add_outgroup(tree, num.shared.var = cloneAge-max(coal_times))
  #tree <- drop.tip(tree, "zeros", collapse.singles = F)

  # Sanity check
  #suppressWarnings(stopifnot(distRoot(tree) == cloneAge))

  # Return the tree created from the coalescence times drawn from Lambert distribution
  return(tree)

}
