#IsingMLE function:

#' IsingMLEmtp2 <to be documented>
#'
#' @param G <to be documented>
#' @param xBar <to be documented>
#' @param M <to be documented>
#' @param data <to be documented>
#' @param epsilon <to be documented>
#' @param maxIter <to be documented>
#' @return <to be documented>
#' @export
IsingMLEmtp2 <- function(G, xBar = NULL, M = NULL, data = NULL, epsilon = 1e-4, maxIter = 100){
  V <- G[[1]]
  E <- G[[2]]

  if (is.null(data)) {
    if (is.null(xBar) | is.null(M)) {
      stop("Must input either a data set or sufficient statistics")
    }
    mu <- xBar
  } else {
    M <- calculateM(as.matrix(data), dim(data)[1], dim(data)[2])
    xBar <- calculatexBar(as.matrix(data), dim(data)[1], dim(data)[2])
    mu <- xBar
  }

  #if(!(all(abs(xBar) < 1))) {
  #  stop(cat("input doesn't fulfill conditions for existance of MLE.\nEmpirical mean contains 1 or -1."))
  #}

  d <- length(mu)

  mCond <- matrix(0, d, d)
  for (i in (1+seq_len(d-1))) {
    for (j in seq_len(i-1))
      mCond[i] <- M[i, j]
  }
  #if(!(all(abs(mCond) < 1))) {
  #  stop(cat("input doesn't fulfill conditions for existance of MLE.\nEmpirical correlation matrix contains 1 or -1."))
  #}

  Xi <- diag(d)

  p <- createP(mu, d)

  ePlus <- createEPlus(E, M, xBar)

  eHat <- matrix(NA, nrow(ePlus), 2)
  eHatOmitNA <- na.omit(eHat)

  condition <- c(Inf)
  condition2 <- calculateCondition2(E, M, Xi)

  empiricalList <- calculateEmpirical(ePlus, M, xBar)
  empirical <- empiricalList$empirical
  containsZeroes <- empiricalList$containsZeroes
  if (ncol(empirical) == 0) {
    stop(cat("input doesn't fulfill conditions for existance of MLE.\n Produced negative values of the empirical distribution"))
  }

  if (containsZeroes | nrow(empirical) == 0) {
    iter <- 0L
    while((iter < maxIter) & (max(abs(mu-xBar)) > epsilon | suppressWarnings(max(condition2) > epsilon) | suppressWarnings(max(condition) > epsilon))){
      pAndEHatResult <- calculateNewPAndEHatBoundary(ePlus, d, empirical, p, eHat)
      p <- pAndEHatResult$p
      eHat <- pAndEHatResult$eHat

      mu <- calculateNewMu(p, d)
      Xi <- calculateNewXi(p, d)

      eHatOmitNA <- na.omit(eHat)

      condition <- calculateCondition(eHatOmitNA, M, Xi);
      condition2 <- calculateCondition2(E, M, Xi);

      iter <- iter + 1L
    }

    if (iter >= maxIter) {
      warning("MLE did not converge; maximum number of iterations reached")
    }
    warning("Maximum likelihood estimate is on the boundary")

    return(list(p_hat = p, G_hat = list(V_hat = V, E_hat = eHatOmitNA), mu_hat = mu, Xi_hat = Xi, number_of_iterations = iter))
  } else {
    iter <- 0L
    capitalJ <- matrix(0, d, d)
    while((iter < maxIter) & (max(abs(mu-xBar)) > epsilon | suppressWarnings(max(condition2) > epsilon) | suppressWarnings(max(condition) > epsilon))){
      pAndEHatResult <- calculateNewPAndEHat(ePlus, d, empirical, p, eHat, capitalJ, iter)
      p <- pAndEHatResult$p
      eHat <- pAndEHatResult$eHat
      capitalJ <- pAndEHatResult$J

      mu <- calculateNewMu(p, d)
      Xi <- calculateNewXi(p, d)

      eHatOmitNA <- na.omit(eHat)

      condition <- calculateCondition(eHatOmitNA, M, Xi)
      condition2 <- calculateCondition2(E, M, Xi)

      iter <- iter + 1L
    }

    h <- calculateH(d, p)

    if (iter >= maxIter) {
      warning("MLE did not converge; maximum number of iterations reached")
    }

    return(list(p_hat = p, G_hat = list(V_hat = V, E_hat = eHatOmitNA), mu_hat = mu, Xi_hat = Xi, h_hat = h, J_hat = capitalJ, number_of_iterations = iter))
  }
}
