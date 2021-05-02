# Function fitting without MTP2 condition:

#' Calculate the maximum likelihood estimator in the Ising model.
#'
#' \code{IsingMLE} fits graphical Ising models using an IPS-type algorithm.
#'
#' Add detailed description.\cr The graph G specifies which undirected graph the maximum likelihood estimator should be Markov with respect to.
#' The graph should be specified as a list with the first element being a vector containing the vertices of the graph and the second element
#' being a matrix containing the edges in the rows. Make sure the naming of the vertices and the edges match.\cr
#' If the maximum likelihood estimator is on the boundary of the parameter space only the mean value parameters will be returned.
#' Otherwise the canonical parameters h and J will also be returned.\cr
#' The index of the probability vector returned by the function corresponds to the binary value of the matching observation with -1 encoded as 0.
#'
#' @param G A list containing a vector with vertices and a matrix containing edges in the rows.
#' @param xBar An optional vector containing the sample first moment of the data. Must be specified if data is unspecified.
#' @param M An optional matrix containing the sample second moment of the data. Must be specified if data is unspecified.
#' @param data An optional matrix or data frame containing observations from d binary variables with outcomes coded as -1 and 1.
#' @param epsilon A numeric value > 0 specifying the tolerance for when the algorithm is considered to have converged.
#' @param maxIter An integer specifying the maximum number of iterations to run the algorithm.
#' @return \code{IsingMLE} returns a list with the estimated distribution, estimated graph, estimated parameters, and number of iterations until the algorithm converged. Expand on this!
#' @export
IsingMLE <- function(G, xBar = NULL, M = NULL, data = NULL, epsilon = 1e-4, maxIter = 100L){
  if (!is.list(G) || length(G) !=2) {stop("G must be a list of length two.")}
  if (!is.vector(G[[1]]) || !is.matrix(G[[2]])) { stop("G must contain a vector with vertices and a matrix with two columns having edges in the rows.")}
  V <- seq_along(G[[1]])
  E <- matrix(data = 0, nrow = nrow(G[[2]]), ncol = ncol(G[[2]]))
  for (i in V) {
    E[G[[2]] == G[[1]][i]] <- i
  }

  if (is.null(data)) {
    if (is.null(xBar) | is.null(M)) {
      stop ("Must input either a data set or sufficient statistics")
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

  #mCond <- matrix(0, d, d)
  #for(i in (1+seq_len(d-1))) {
  #  for(j in seq_len(i-1))
  #    mCond[i] <- M[i, j]
  #}
  #if(!(all(abs(mCond) < 1))) {
  #  stop(cat("input doesn't fulfill conditions for existance of MLE.\nEmpirical correlation matrix contains 1 or -1."))
  #}

  Xi <- diag(d)

  p <- createP(mu, d)

  eHat <- matrix(NA, nrow(E), 2)
  eHatOmitNA <- na.omit(eHat)

  condition <- c(Inf)

  empiricalList <- calculateEmpirical(E, M, xBar)
  empirical <- empiricalList$empirical
  containsZeroes <- empiricalList$containsZeroes
  if (ncol(empirical) == 0) {
    stop (cat("input doesn't fulfill conditions for existance of MLE.\n Produced negative values of the empirical distribution"))
  }

  if (containsZeroes) {
    iter <- 0L
    while ((iter < maxIter) & (max(abs(mu-xBar)) > epsilon | max(condition) > epsilon)){
      pAndEHatResult <- calculateNewPAndEHatnoMTP2Boundary(E, d, empirical, p, eHat)
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
      warning ("MLE did not converge; maximum number of iterations reached")
    }
    warning ("Maximum likelihood estimate is on the boundary")

    for (i in V) {
      eHatOmitNA[eHatOmitNA == V[i]] <- G[[1]][i]
    }
    attr(eHatOmitNA, "na.action") <- NULL
    attr(eHatOmitNA, "class") <- NULL

    return (list(p_hat = p, G_hat = list(V_hat = G[[1]], E_hat = eHatOmitNA), mu_hat = mu, Xi_hat = Xi, number_of_iterations = iter))
  } else {
    iter <- 0L
    capitalJ <- matrix(0, d, d)
    while ((iter < maxIter) & (max(abs(mu-xBar)) > epsilon | max(condition) > epsilon)){
      pAndEHatResult <- calculateNewPAndEHatnoMTP2(E, d, empirical, p, eHat, capitalJ, iter)
      p <- pAndEHatResult$p
      eHat <- pAndEHatResult$eHat
      capitalJ <- pAndEHatResult$J

      mu <- calculateNewMu(p, d)
      Xi <- calculateNewXi(p, d)

      eHatOmitNA <- na.omit(eHat)

      condition <- calculateCondition(eHatOmitNA, M, Xi);

      iter <- iter + 1L
    }

    h <- calculateH(d, p);

    if (iter >= maxIter) {
      warning ("MLE did not converge; maximum number of iterations reached")
    }

    for (i in V) {
      eHatOmitNA[eHatOmitNA == V[i]] <- G[[1]][i]
    }
    attr(eHatOmitNA, "na.action") <- NULL
    attr(eHatOmitNA, "class") <- NULL

    return (list(p_hat = p, G_hat = list(V_hat = G[[1]], E_hat = eHatOmitNA), mu_hat = mu, Xi_hat = Xi, h_hat = h, J_hat = capitalJ, number_of_iterations = iter))
  }
}
