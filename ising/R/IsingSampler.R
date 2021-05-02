library(Rcpp)

#' Convert integer to binary with 0 represented with -1 and 1 represented with 1
#'
#' \code{IntToSign} converts an integer to a binary number where 0 is represented with -1 and 1 is represented with 1.
#'
#' @param x Integer
#' @param digits Integer specifying how many digits should be returned
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

#' IsingSampler <to be documented>
#'
#' @param h <to be documented>
#' @param J <to be documented>
#' @param p <to be documented>
#' @param N <to be documented>
#' @param obs <to be documented>
#' @param int <to be documented>
#' @param matrix <to be documented>
#' @return <to be documented>
#' @export
IsingSampler <- function(h = NULL, J = NULL, p = NULL, N = 1000, obs = TRUE, int = FALSE, matrix = FALSE) {
  if(!is.null(h)) {
    if(!(is.numeric(h))) stop("h must be a numeric vector")
  }
  if(!is.null(J)) {
    if(!(is.matrix(J) & isSymmetric(J) & length(h) == dim(J)[1])) stop("J must be a symmetric matrix with dimensions equal to the lenght of h.")
  }

  if(is.null(p)) {
    d <- length(h)
    p <- computeP(h, J)
  } else {
    d <- log2(length(p))
  }

  Sample <- matrix(0, N, d)

  IntSample <- sample(x = 0:(2^d-1), size = N, replace = TRUE, prob = p)

  for(i in 1:N) Sample[i, ] <- IntToSign(IntSample[i], d)

  if (matrix) {
    if(int & obs) {
      return(list(observations = Sample, observations_int = IntSample))
    } else if(obs) {
      return(Sample)
    } else if(int) {
      return(IntSample)
    } else {
      stop("At least one of obs or int should be set to TRUE")
    }
  } else {
    if(int & obs) {
      return(list(observations = as.data.frame(Sample), observations_int = IntSample))
    } else if(obs) {
      return(as.data.frame(Sample))
    } else if(int) {
      return(IntSample)
    } else {
      stop("At least one of obs or int should be set to TRUE")
    }
  }
}
