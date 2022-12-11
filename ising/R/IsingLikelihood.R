# IsingLogLikelihood function:

#' Calculate log-likelihood in Ising models
#'
#' \code{IsingLogLikelihood} calculates the log-likelihood of a data set for a
#' specified Ising model.
#'
#' The function must be provided with either the natural parameters h and J or
#' the distribution p. \cr
#' A probability in the probability p belongs to the observation whose binary
#' value matches the value of the index (using zero-indexing).
#'
#' @param data Matrix or data frame containing samples from d binary variables.
#' @param h Vector containing the canonical parameters h.
#' @param J Matrix containing the canonical parameters J.
#' @param p Vector containing the distribution of the Ising model.
#' @return \code{IsingLogLikelihood} returns a numeric value with the
#' log-likelihood of the data set for the specified Ising model.
#' @export
IsingLogLikelihood <- function(data,
                               h = NULL,
                               J = NULL,
                               p = NULL) {
  # the functions call lower level c++ functions based on whether h and J or p is specified.
  if(is.null(h) & is.null(J) & !is.null(p)) {
    return(isingLogLikelihoodProbs(data, p))
  }
  if(!is.null(h) & !is.null(J) & is.null(p)) {
    return(IsingLogLikelihood <- isingLogLikelihoodhJ(data, h, J))
  }
  stop("Must specify either the natural parameters h and J or the probability vector p")
}
