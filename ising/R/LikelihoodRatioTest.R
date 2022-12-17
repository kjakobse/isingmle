# likelihoodRatioTest function:

#' Calculate values of the log-likelihood ratio test statistic using simulated
#' data.
#'
#' \code{likelihoodRatioTest} simulates data sets from a specified
#'  null-distribution and for each calculates the log-likelihood ratio test
#'  statistic for the test against the assumed model.
#'
#' The function uses the parallel implementation in \code{parSapply} to allow
#' running multiple simulations at a time. The parameter nCores should not
#' exceed the number of logical cores on the host system. For larger values of
#' d available memory also needs to be considered since the memory requirements
#' scale linearly with the value of nCores. \cr As is standard for hypothesis
#' testing the large model should contain the null model.
#'
#' @param nSim Integer specifying number of simulations
#' @param p Vector containing the distribution of the null model.
#' @param d Integer specifying number of binary variables.
#' @param N Integer specifying number of observations to be simulated from the
#' null model.
#' @param largeModel Character string specifying the assumed model. Can be
#' "IsingMTP2", "Ising", or "FullBinary".
#' @param nullModel Character string specifying the null model. Can be
#' "IsingMTP2" or "Ising". The null model should always be contained in
#' largeModel.
#' @param epsilon A numeric value > 0 specifying the tolerance for when the
#' fitting algorithm is considered to have converged.
#' @param GLarge List containing vector of vertices and matrix of edges.
#' Must be specified when largeModel is "IsingMTP2" or "Ising".
#' @param GNull List containing vector of vertices and matrix of edges.
#' @param maxIter Integer specifying maximum number of iterations to run the
#' fitting algorithm.
#' @param nCores Integer specifying the number of cores to use for the
#' simulation.
#' @param zeroReplace A boolean value indicating whether or not to replace
#' zeroes in the empirical distribution.
#' @param ReplaceValue A numeric value >0 specifying which value to replace
#' zeroes with.
#' @return \code{likelihoodRatioTest} returns a vector of length nSim with the
#' likelihood ratio test statistics for the data set simulated under the null
#' model.
#' @export
likelihoodRatioTest <- function(nSim,
                                p,
                                d,
                                N,
                                largeModel,
                                nullModel,
                                epsilon,
                                GLarge = NULL,
                                GNull = NULL,
                                maxIter = 100L,
                                nCores = 1L,
                                zeroReplace = FALSE,
                                ReplaceValue = 1e-10) {
  if (!is.numeric(nSim) || ((nSim %% 1) != 0) || (nSim < 1)) {
    stop("nSim must be an integer of value greater than zero")
  }
  # If nCores is one, bootstrapLikelihoodRatio is called with replicate to
  # create a vector of log-likelihood ratio test statistics.
  # If nCores is larger than one, bootstrapLikelihoodRatio is called with
  # parSapply to create a vector of log-likelihood ratio test statistics.
  if (nCores == 1) {
    res <- replicate(
      nSim,
      bootstrapLikelihoodRatio(
        p = p,
        d = d,
        N = N,
        largeModel = largeModel,
        nullModel = nullModel,
        epsilon = epsilon,
        GLarge = GLarge,
        GNull = GNull,
        maxIter = maxIter,
        zeroReplace = zeroReplace,
        ReplaceValue = ReplaceValue
      )
    )
  } else {
    cl <- makeCluster(nCores)
    res <- parSapply(
      cl,
      integer(nSim),
      function(...) {
        bootstrapLikelihoodRatio(
          p = p,
          d = d,
          N = N,
          largeModel = largeModel,
          nullModel = nullModel,
          epsilon = epsilon,
          GLarge = GLarge,
          GNull = GNull,
          maxIter = maxIter,
          zeroReplace = zeroReplace,
          ReplaceValue = ReplaceValue
        )
      },
      p = p,
      d = d,
      N = N,
      largeModel = largeModel,
      nullModel = nullModel,
      epsilon = epsilon,
      GLarge = GLarge,
      GNull = GNull,
      maxIter = maxIter,
      zeroReplace = zeroReplace,
      ReplaceValue = ReplaceValue
    )
    stopCluster(cl)
  }

  return(res)
}
