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
#' @importFrom ape "is.ultrametric"
#' @examples
#' data(exampleTrees)
#' df <- moments(exampleTrees[[1]])
moments <- function(subtree, alpha = 0.05){

  ptm <- proc.time()

  # Must be of class phylo
  if(class(subtree) != "phylo"){
    stop("Tree must be of class phylo. Use as.phylo function to convert if the
    formatting is correct. Otherwise, see ape documentation
    https://cran.r-project.org/web/packages/ape/ape.pdf")
  }

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

  moments_growth_rate <- (pi/sqrt(3))*1/(sd(branching.times(t1)))
  moments_growth_rate_lb <- moments_growth_rate*sqrt(1+4*qnorm(alpha/2)/sqrt(5*n))
  moments_growth_rate_ub <- moments_growth_rate*sqrt(1-4*qnorm(alpha/2)/sqrt(5*n))

  moments <- c("lb_2.5" = moments_growth_rate_lb, "estimate" = moments_growth_rate,
               "ub_97.5" = moments_growth_rate_ub)

  runtime <- proc.time() - ptm

  return(data.frame("lowerBound" = growthRate_lb, "estimate" = growthRate,
                    "upperBound" = growthRate_ub, "sumInternalLengths" = IL,
                    "n" = n, "alpha" = alpha, "hasStem" = hasStem,
                    "includeStem" = includeStem, "runtime_s" = runtime))
}
