# Function fitting without MTP2 condition:

#' Calculate the maximum likelihood estimator in the Ising model without MTP2
#' constraints.
#'
#' \code{IsingMLE} fits graphical Ising models using an IPS-type algorithm.
#'
#' The graph G specifies which undirected graph the maximum likelihood estimator
#' should be Markov with respect to.
#' The graph should be specified as a list with the first element being a
#' vector containing the vertices of the graph and the second element
#' being a matrix containing the edges in the rows. Make sure the naming of
#' the vertices and the edges match.\cr
#' If the maximum likelihood estimator is on the boundary of the parameter space
#' only the mean value parameters will be returned.
#' Otherwise the canonical parameters h and J will also be returned.\cr
#' A probability in the probability vector returned by the function belongs to
#' the observation whose binary value matches the value of the index
#' (using zero-indexing).
#'
#' @param G A list containing a vector with vertices and a matrix containing
#' edges in the rows.
#' @param xBar An optional vector containing the sample first moment of the
#' data. Must be specified if data is unspecified.
#' @param M An optional matrix containing the sample second moment of the data.
#' Must be specified if data is unspecified.
#' @param data An optional matrix or data frame containing samples from d binary
#' variables with outcomes coded as -1 and 1.
#' @param epsilon A numeric value >0 specifying the tolerance for when the
#' algorithm is considered to have converged.
#' @param maxIter An integer specifying the maximum number of iterations to run
#' the algorithm.
#' @param zeroReplace A boolean value indicating whether or not to replace
#' zeroes in the empirical distribution.
#' @param ReplaceValue A numeric value >0 specifying which value to replace
#' zeroes with.
#' @return \code{IsingMLE} returns a list with the estimated distribution,
#' estimated graph, estimated parameters, and number of iterations until the
#' algorithm converged.
#' @author Kim Daniel Jakobsen
#' @examples
#' 1 + 1
#'
#' @export
#' @import stats
IsingMLE <- function(
    G,
    xBar = NULL,
    M = NULL,
    data = NULL,
    epsilon = 1e-4,
    maxIter = 100L,
    zeroReplace = FALSE,
    ReplaceValue = 1e-10) {
  # Encode the vertices in G as the integers from 1 to d:
  if (!is.list(G) || length(G) != 2) {
    stop("G must be a list of length two.")
  }
  if (!is.vector(G[[1]]) || !is.matrix(G[[2]])) {
    stop(
      "G must contain a vector with vertices and a matrix ",
      "with two columns having edges in the rows."
    )
  }
  V <- seq_along(G[[1]])
  E <- matrix(
    data = 0,
    nrow = nrow(G[[2]]),
    ncol = ncol(G[[2]])
  )
  for (i in V) {
    E[G[[2]] == G[[1]][i]] <- i
  }

  # Initiate the mean value parameter mu with the empirical value. If not given
  # directly the empirical moments are calculated from the provided data set:
  if (!is.null(data) && (is.null(xBar) || is.null(M))) {
    M <- calculateM(as.matrix(data), dim(data)[1], dim(data)[2])
    xBar <- calculatexBar(as.matrix(data), dim(data)[1], dim(data)[2])
    mu <- xBar
  }
  if (!is.null(xBar) && !is.null(M)) {
    mu <- xBar
  }
  if (!exists("mu", mode = "double")) {
    stop("Must input either a data set or sufficient statistics")
  }

  # initiate d, the number of binary variables:
  d <- length(mu)

  # initiate the mean value parameter Xi:
  Xi <- diag(d)

  # initiate the distribution given by mu and Xi:
  p <- createP(mu, d)

  # check if the initial distribution contains zeroes:
  if (!(all(p > 1e-15))) {
    pcontainszeroes <- 1
  } else {
    pcontainszeroes <- 0
  }

  # initiate a matrix for the edges of the fitted graph G_hat:
  eHat <- matrix(NA, nrow(E), 2)
  eHatOmitNA <- na.omit(eHat)

  # Initiate the condition on Xi. Inf ensures that the statement
  # max(condition > epsilon) is always true initially:
  # condition checks if difference between Xi and M for entries in E are zero.
  condition <- Inf

  if (zeroReplace) {
    # Calculate the empirical distribution for each variable pair in E and
    # replace zeroes with ReplaceValue
    empirical <- calculateEmpiricalReplaceZeroes(
      E,
      M,
      xBar,
      ReplaceValue
    )
    econtainsZeroes <- 0
    if (ncol(empirical) == 0) {
      stop(cat(
        "input doesn't fulfill conditions for existance of MLE.",
        "\n",
        "Produced negative values of the empirical distribution",
        sep = ""
      ))
    }
  } else {
    # Calculate the empirical distribution for each variable pair in E and
    # check if any contain zeroes:
    empiricalList <- calculateEmpirical(E, M, xBar)
    empirical <- empiricalList$empirical
    econtainsZeroes <- empiricalList$containsZeroes
    if (ncol(empirical) == 0) {
      stop(cat(
        "input doesn't fulfill conditions for existance of MLE.",
        "\n",
        "Produced negative values of the empirical distribution.",
        sep = ""
      ))
    }
  }

  # If the empirical distribution contains zeroes a version of the algorithm
  # converging on the boundary is used:
  if (econtainsZeroes || nrow(empirical) == 0 || pcontainszeroes) {
    iter <- 0L
    # While loop updating the distribution until the convergence criteria is
    # meet or the max number of iterations is reached:
    while (
      (iter < maxIter) &&
        (
          max(abs(mu - xBar)) > epsilon ||
            max(condition) > epsilon
        )
    ) {
      # Update the distribution for each variable pair in E, and update the
      # edges in the fitted graph:
      pAndEHatResult <- calculateNewPAndEHatnoMTP2Boundary(
        E,
        d,
        empirical,
        p
      )
      p <- pAndEHatResult$p

      # Calculate the updated values of the mean value parameters:
      mu <- calculateNewMu(p, d)
      Xi <- calculateNewXi(p, d)

      # Calculate the updated condition in the convergence criteria:
      condition <- calculateCondition(E, M, Xi)

      iter <- iter + 1L
    }

    if (iter >= maxIter) {
      warning("MLE did not converge; maximum number of iterations reached")
    }
    warning("Maximum likelihood estimate is on the boundary")

    # Return the fitted distribution, graph, mean value parameters and number
    # of iterations until convergence:
    return(
      list(
        p_hat = p,
        mu_hat = mu,
        Xi_hat = Xi,
        number_of_iterations = iter
      )
    )
  } else {
    iter <- 0L
    # Initiate a matrix for the canonical parameter J:
    capitalJ <- matrix(0, d, d)
    # While loop updating the distribution until the convergence criteria is
    # meet or the max number of iterations is reached:
    while (
      (iter < maxIter) &&
        (
          max(abs(mu - xBar)) > epsilon ||
            max(condition) > epsilon
        )
    ) {
      # Update the distribution for each variable pair in E, update the edges
      # in the fitted graph, and calculate the value of the canonical parameter
      # J for the updated Ising model:
      pAndEHatResult <- calculateNewPAndEHatnoMTP2(
        E,
        d,
        empirical,
        p,
        eHat,
        capitalJ
      )
      p <- pAndEHatResult$p
      eHat <- pAndEHatResult$eHat
      capitalJ <- pAndEHatResult$J

      # Calculate the updated values of the mean value parameters:
      mu <- calculateNewMu(p, d)
      Xi <- calculateNewXi(p, d)

      eHatOmitNA <- na.omit(eHat)

      # Calculate the updated condition in the convergence criteria:
      condition <- calculateCondition(E, M, Xi)

      iter <- iter + 1L
    }

    # Calculate the canonical parameter h for the fitted Ising model:
    h <- calculateH(d, p)

    if (iter >= maxIter) {
      warning("MLE did not converge; maximum number of iterations reached")
    }

    # Recode the vertices in G to the names specified in the input:
    for (i in V) {
      eHatOmitNA[eHatOmitNA == V[i]] <- G[[1]][i]
    }
    attr(eHatOmitNA, "na.action") <- NULL
    attr(eHatOmitNA, "class") <- NULL

    # return the fitted distribution, graph, mean value parameters,
    # the canonical parameters, and number of iterations until convergence:
    return(
      list(
        p_hat = p,
        G_hat = list(V_hat = G[[1]], E_hat = eHatOmitNA),
        mu_hat = mu,
        Xi_hat = Xi,
        h_hat = h,
        J_hat = capitalJ,
        number_of_iterations = iter
      )
    )
  }
}
