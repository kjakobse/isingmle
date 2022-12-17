# IntToSign function:

#' Internal function. Convert integer to binary with 0 represented with -1 and 1
#' represented with 1.
#'
#' \code{IntToSign} converts an integer to a binary number where 0 is
#' represented with -1 and 1 is represented with 1.
#'
#' @param x Integer to be converted.
#' @param digits Integer specifying how many digits should be returned.
#' @return \code{IntToSign} returns a vector with length equal to digits,
#' containing the sign representation of the input integer.
IntToSign <- function(x, digits) {
  i <- 0L
  Sign <- rep(-1, digits)
  while (x > 0) {
    if(x %% 2L == 1) {
      Sign[digits - i] <- 1
    }
    x <- x %/% 2L
    i <- i + 1L
  }
  Sign
}
