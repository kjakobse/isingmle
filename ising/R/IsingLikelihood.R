#' IsingLogLikelihood <to be documented>
#'
#' @param data <to be documented>
#' @param h <to be documented>
#' @param J <to be documented>
#' @param p <to be documented>
#' @return <to be documented>
#' @export
IsingLogLikelihood <- function(data, h = NULL, J = NULL, p = NULL) {
  if(is.null(h) & is.null(J) & !is.null(p)) {
    IsingLogLikelihood <- isingLogLikelihoodProbs(data, p)
  } else if(!is.null(h) & !is.null(J) & is.null(p)) {
    IsingLogLikelihood <- isingLogLikelihoodhJ(data, h, J)
  } else {
    stop("Must specify either the natural parameters h and J or the probability vector p")
  }
  IsingLogLikelihood
}
