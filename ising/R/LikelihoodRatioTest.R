library(parallel)

#' Calculate values of the likelihood ratio test statistic using simulated data.
#'
#' \code{LikelihoodRatioTest} simulates data sets from a specified null-distribution and for each calculates the likelihood ratio test statistic for the test against the assumed model.
#'
#' Detailed description.
#'
#' @param nSim Integer specifying number of simulations
#' @param p Vector containing the distribution of the null model.
#' @param d Integer specifying number of binary variables.
#' @param N Integer specifying number of observations to be simulated from the null model.
#' @param largeModel Character string specifying the assumed model. Can be "IsingMTP2", "Ising", or "FullBinary".
#' @param nullModel Character string specifying the null model. Can be "IsingMTP2" or "Ising". The null model should always be contained in largeModel.
#' @param epsilon A numeric value > 0 specifying the tolerance for when the fitting algorithm is considered to have converged.
#' @param GLarge List containing vector of vertices and matrix of edges. Must be specified when largeModel is "IsingMTP2" or "Ising".
#' @param GNull List containing vector of vertices and matrix of edges.
#' @param maxIter Integer specifying maximum number of iterations to run the fitting algorithm.
#' @param nCores Integer specifying the number of cores to use for the simulation.
#' @return \code{LikelihoodRatioTest} returns a vector of length nSim with the likelihood ratio test statistics for the data set simulated under the null model.
#' @export
likelihoodRatioTest <- function(nSim, p, d, N, largeModel, nullModel, epsilon, GLarge = NULL, GNull = NULL, maxIter = 100L, nCores = 1L) {
  if (!is.numeric(nSim) || ((nSim %% 1) != 0) || (nSim < 1)) {
    stop("nSim must be an integer of value greater than zero")
  }
  nList <- sapply(1:nSim, list)
  if (nCores == 1) {
    res <- sapply(nList, bootstrapLikelihoodRatio, p = p, d = d, N = N, largeModel = largeModel, nullModel = nullModel, epsilon = epsilon, GLarge = GLarge, GNull = GNull, maxIter = maxIter)
  } else {
    cl <- makeCluster(nCores)
    res <- parSapply(cl, nList, bootstrapLikelihoodRatio, p = p, d = d, N = N, largeModel = largeModel, nullModel = nullModel, epsilon = epsilon, GLarge = GLarge, GNull = GNull, maxIter = maxIter)
    stopCluster(cl)
  }

  return(res)
}
