library(Rcpp)

# IntToSign function:

#' Convert integer to binary with 0 represented with -1 and 1 represented with 1
#'
#' \code{IntToSign} converts an integer to a binary number where 0 is represented with -1 and 1 is represented with 1.
#'
#' @param x Integer to be converted.
#' @param digits Integer specifying how many digits should be returned.
#' @return \code{IntToSign} returns a vector with length equal to digits, containing the sign representation of the input integer.
#' @export
IntToSign <- function(x, digits) {
  i <- 0L
  Sign <- rep(-1, digits)
  while (x > 0) {
    if(x %% 2L == 0) {
      Sign[digits - i] <- -1
    }
    else {
      Sign[digits - i] <- 1
    }
    x <- x %/% 2L
    i <- i + 1L
  }
  Sign
}

# IsingSampler function:

#' Generate data sets from the Ising model
#'
#' \code{IsingSampler} generates a data set of a specified length from an Ising model.
#'
#' The matrix J should be symmetric and have zeroes in the diagonal. \cr
#' A probability in the probability vector p is assumed to belong to the observation whose binary value matches the value of the index (using zero-indexing).
#'
#' @param h Vector, specifying the canonical parameter h from the Ising model. Optional if p is specified.
#' @param J Matrix, specifying the canonical parameter J from the Ising model. Optional if p is specified.
#' @param p Vector, specifying the distribution of the Ising model. Optional if h and J are specified.
#' @param N Integer, specifying the number of observations in the generated data set.
#' @param obs Logical, if TRUE the observations are returned in binary form.
#' @param int Logical, if TRUE, the observations are returned as integers.
#' @param matrix Logical, if TRUE the binary form of the observations are returned in a matrix. If FALSE the observations are returned in a data frame.
#' @return \code{IsingSampler} returns a simulated data set of length N from the Ising model specified by h and J or by p. The observations are returned in binary form if the obs parameter is TRUE, and as integers if the parameter int is TRUE. The matrix parameter controls if the observations in binary form are returned in a matrix or in a data frame.
#' @export
IsingSampler <- function(h = NULL, J = NULL, p = NULL, N = 1000L, obs = TRUE, int = FALSE, matrix = FALSE) {
  if(!is.null(h)) {
    if(!(is.numeric(h))) stop("h must be a numeric vector")
  }
  if(!is.null(J)) {
    if(!(is.matrix(J) & isSymmetric(J) & length(h) == dim(J)[1])) stop("J must be a symmetric matrix with dimensions equal to the lenght of h.")
  }

  # Initiate the number of binary variables d, and the distribution p:
  if(is.null(p)) {
    d <- length(h)
    p <- computeP(h, J)
  } else {
    d <- log2(length(p))
  }

  # Sample indices whose binary representation corresponds to the observations in the data set:
  IntSample <- sample(x = 0:(2^d-1), size = N, replace = TRUE, prob = p)

  # Return the sampled data in the format specified by the input. Either just the integers, just the corresponding binary observations, or both:
  if (matrix) {
    if(int & obs) {
      Sample <- matrix(0, N, d)
      for(i in 1:N) Sample[i, ] <- IntToSign(IntSample[i], d)
      return(list(observations = Sample, observations_int = IntSample))
    } else if(obs) {
      Sample <- matrix(0, N, d)
      for(i in 1:N) Sample[i, ] <- IntToSign(IntSample[i], d)
      return(Sample)
    } else if(int) {
      return(IntSample)
    } else {
      stop("At least one of obs or int should be set to TRUE")
    }
  } else {
    if(int & obs) {
      Sample <- matrix(0, N, d)
      for(i in 1:N) Sample[i, ] <- IntToSign(IntSample[i], d)
      return(list(observations = as.data.frame(Sample), observations_int = IntSample))
    } else if(obs) {
      Sample <- matrix(0, N, d)
      for(i in 1:N) Sample[i, ] <- IntToSign(IntSample[i], d)
      return(as.data.frame(Sample))
    } else if(int) {
      return(IntSample)
    } else {
      stop("At least one of obs or int should be set to TRUE")
    }
  }
}
