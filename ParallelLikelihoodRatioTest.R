library(ising)

# p, d, N, largeModel, nullModel, epsilon, GLarge = NULL, GNull = NULL, maxiter = 1000
nSim <- 10000
p <- rep(1/16, 16)
d <- 4
N <- 1000
largeModel <- "Ising"
nullModel <- "IsingMtp2"
epsilon <- 1e-4
GLarge <- list(1:4, matrix(c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4), 6, 2, T))
GNull <- list(1:4, matrix(c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4), 6, 2, T))
maxIter <- 1000
nCores <- 6

likelihoodRatioTest(nSim = nSim, p = p, d = d, N = N, largeModel = largeModel, nullModel = nullModel, epsilon = epsilon, GLarge = GLarge, GNull = GNull, maxIter = maxIter, nCores = nCores)
