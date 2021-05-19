# bootstrapLikelihoodRatio function:

#' Simulate value of likelihood ratio test statistic
#'
#' \code{bootstrapLikelihoodRatio} calculates the likelihood ratio test statistic from a data set simulated from the specified null-distribution.
#'
#' This function is intended to be called by \code{LikelihoodRatioTest} and not called directly.
#'
#' @param simNum Integer
#' @param p Vector containing the distribution of the null model.
#' @param d Integer specifying number of binary variables.
#' @param N Integer specifying number of observations to be simulated from the null model.
#' @param largeModel Character vector specifying the assumed model. Can be "IsingMTP2", "Ising", or "FullBinary".
#' @param nullModel Character vector specifying the null model. Can be "IsingMTP2" or "Ising". The null model should always be contained in largeModel.
#' @param epsilon A numeric value > 0 specifying the tolerance for when the fitting algorithm is considered to have converged.
#' @param GLarge List containing vector of vertices and matrix of edges. Must be specified when largeModel is "IsingMTP2" or "Ising".
#' @param GNull List containing vector of vertices and matrix of edges.
#' @param maxIter integer specifying maximum number of iterations to run the fitting algorithm.
#' @param zeroReplace A boolean value indicating whether or not to replace zeroes in the empirical distribution.
#' @param ReplaceValue A numeric value >0 specifying which value to replace zeroes with.
#' @return \code{bootstrapLikelihoodRatio} returns a numeric value with the likelihood ratio test statistic for the data set simulated under the null model.
#' @export
bootstrapLikelihoodRatio <- function(simNum, p, d, N, largeModel, nullModel, epsilon, GLarge = NULL, GNull = NULL, maxIter = 100L, zeroReplace = FALSE, ReplaceValue = 1e-10) {
  # Generate a sample from the given null model with N observations:
  Listsample <- IsingSampler(p = p, N = N, obs = TRUE, int = TRUE, matrix = TRUE)
  sample <- Listsample$observations
  intSample <- Listsample$observations_int
  # From the simulated data, calculate the log-likelihood for the specified null model:
  if (nullModel == "IsingMTP2") {
    if (is.null(GNull)) {
      stop ("GNull must be specified when nullModel is IsingMTP2")
    }
    fit_null <- IsingMLEmtp2(G = GNull, data = sample, epsilon = epsilon, maxIter = maxIter, zeroReplace = zeroReplace, ReplaceValue = ReplaceValue)
    LogLikelihoodNull <- IsingLogLikelihood(data = sample, p = fit_null$p_hat)
  } else if (nullModel == "Ising") {
    if (is.null(GNull)) {
      stop ("GNull must be specified when nullModel is Ising")
    }
    fit_null <- IsingMLE(G = GNull, data = sample, epsilon = epsilon, maxIter = maxIter, zeroReplace = zeroReplace, ReplaceValue = ReplaceValue)
    LogLikelihoodNull <- IsingLogLikelihood(data = sample, p = fit_null$p_hat)
  } else {
    stop ("nullModel must be [\"Ising\" or \"IsingMTP2\"]")
  }
  # From the simulated data, calculate the log-likelihood for the specified large model:
  if (largeModel == "Ising") {
    if (is.null(GLarge)) {
      stop ("GLarge must be specified when largeModel is Ising")
    }
    fit_large <- IsingMLE(G = GLarge, data = sample, epsilon = epsilon, maxIter = maxIter, zeroReplace = zeroReplace, ReplaceValue = ReplaceValue)
    LogLikelihoodLarge <- IsingLogLikelihood(data = sample, p = fit_large$p_hat)
  } else if (largeModel == "IsingMTP2") {
    if (is.null(GLarge)) {
      stop ("GLarge must be specified when nullModel is IsingMTP2")
    }
    fit_large <- IsingMLEmtp2(G = GLarge, data = sample, epsilon = epsilon, maxIter = maxIter, zeroReplace = zeroReplace, ReplaceValue = ReplaceValue)
    LogLikelihoodLarge <- IsingLogLikelihood(data = sample, p = fit_large$p_hat)
  } else if (largeModel == "FullBinary") {
    LogLikelihoodLarge <- binaryLogLikelihood(intData = intSample, d = d)
  } else {
    stop ("largeModel must be [\"Ising\", \"IsingMTP2\", or \"FullBinary\"]")
  }

  # return the value of the log-likelihood ratio test statistic:
  return (-2 * (LogLikelihoodNull - LogLikelihoodLarge))
}
