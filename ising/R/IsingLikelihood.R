#' Calculate log-likelihood in Ising models
#'
#' \code{IsingLogLikelihood} calculates the log-likelihood of a data set for a specified Ising model.
#'
#' Detailed description. The function must be provided with either the natural parameters h and J or the distribution p.
#' the i'th entry of p must contain the probability of the observation which has binary value equal to i when -1 is encoded as 0.
#'
#' @param data Matrix or data frame containing samples from d binary variables.
#' @param h Vector containing the canonical parameters h.
#' @param J Matrix containing the canonical parameters J.
#' @param p Vector containing the distribution of the Ising model.
#' @return \code{IsingLogLikelihood} returns a numeric value with the log-likelihood of the data set for the specified Ising model.
#' @export
IsingLogLikelihood <- function(data, h = NULL, J = NULL, p = NULL) {
  # the functions calls lower level c++ functions based on whether h and J or p is specified.
  if(is.null(h) & is.null(J) & !is.null(p)) {
    IsingLogLikelihood <- isingLogLikelihoodProbs(data, p)
  } else if(!is.null(h) & !is.null(J) & is.null(p)) {
    IsingLogLikelihood <- isingLogLikelihoodhJ(data, h, J)
  } else {
    stop("Must specify either the natural parameters h and J or the probability vector p")
  }
  IsingLogLikelihood
}
