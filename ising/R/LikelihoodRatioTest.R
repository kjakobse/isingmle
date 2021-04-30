library(parallel)

#' likelihoodRatioTest <to be documented>
#'
#' @param nSim <to be documented>
#' @param p <to be documented>
#' @param d <to be documented>
#' @param N <to be documented>
#' @param largeModel <to be documented>
#' @param nullModel <to be documented>
#' @param epsilon <to be documented>
#' @param GLarge <to be documented>
#' @param GNull <to be documented>
#' @param maxIter <to be documented>
#' @param nCores <to be documented>
#' @return <to be documented>
#' @export
likelihoodRatioTest <- function(nSim, p, d, N, largeModel, nullModel, epsilon, GLarge = NULL, GNull = NULL, maxIter = 1000, nCores = 1) {
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
