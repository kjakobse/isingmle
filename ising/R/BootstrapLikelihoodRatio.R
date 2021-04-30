#' bootstrapLikelihoodRatio <to be documented>
#'
#' @param simNum <to be documented>
#' @param p <to be documented>
#' @param d <to be documented>
#' @param N <to be documented>
#' @param largeModel <to be documented>
#' @param nullModel <to be documented>
#' @param epsilon <to be documented>
#' @param GLarge <to be documented>
#' @param GNull <to be documented>
#' @param maxIter <to be documented>
#' @return <to be documented>
#' @export
bootstrapLikelihoodRatio <- function(simNum, p, d, N, largeModel, nullModel, epsilon, GLarge = NULL, GNull = NULL, maxIter = 1000) {
  Listsample <- IsingSampler(p = p, N = N, obs = TRUE, int = TRUE, matrix = TRUE)
  sample <- Listsample$observations
  intSample <- Listsample$observations_int
  if (nullModel == "IsingMtp2") {
    if (is.null(GNull)) {
      stop("GNull must be specified when nullModel is IsingMtp2")
    }
    fit_mtp2 <- IsingMLEmtp2(G = GNull, data = sample, epsilon = epsilon, maxIter = maxIter)
    LogLikelihoodNull <- IsingLogLikelihood(data = sample, p = fit_mtp2$p_hat)
  } else {
    stop("nullModel must be [\"IsingMtp2\"]")
  }
  if (largeModel == "Ising") {
    if (is.null(GNull)) {
      stop("GLarge must be specified when largeModel is Ising")
    }
    fit_nomtp2 <- IsingMLE(G = GLarge, data = sample, epsilon = epsilon, maxIter = maxIter)
    LogLikelihoodLarge <-IsingLogLikelihood(data = sample, p = fit_nomtp2$p_hat)
  } else if (largeModel == "FullBinary") {
    LogLikelihoodLarge <- binaryLogLikelihood(intData = intSample, d = d)
  } else {
    stop("largeModel must be [\"Ising\", \"FullBinary\"]")
  }
  cat("simNum = ", simNum, "\n")

  return(-2 * (LogLikelihoodNull - LogLikelihoodLarge))
}
