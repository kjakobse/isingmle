library(parallel)

setwd("D:/development/isingmle")
source("IsingMLECommon.R")
source("IsingMLE.R")
source("IsingMLEmtp2.R")
source("IsingLikelihood.R")
source("BinaryLikelihood.R")
source("IsingSampler.R")
source("BootstrapLikelihoodRatio.R")

# p, d, N, largeModel, nullModel, epsilon, GLarge = NULL, GNull = NULL, maxiter = 1000
p <- rep(1/16, 16)
d <- 4
N <- 1000
largeModel <- "Ising"
nullModel <- "IsingMtp2"
epsilon <- 1e-4
GLarge <- list(1:4, matrix(c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4), 6, 2, T))
GNull <- list(1:4, matrix(c(1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4), 6, 2, T))

cl <- makeCluster(4)
clusterEvalQ(cl, {
  setwd("D:/development/isingmle")
  source("BinaryLikelihood.R")
})
clusterExport(cl, c("IsingSampler", "IntToSign", "computeP", "IsingMLE", "calculateNewPAndEHatnoMTP2", "calculateNewPAndEHatnoMTP2Boundary", "IsingMLEmtp2", "calculateNewPAndEHat", "calculateNewPAndEHatBoundary", "calculatexBar", "calculateM", "createP", "calculateNewMu", "calculateNewXi", "createEPlus", "calculateCondition", "calculateCondition2", "CalculateEmpirical", "calculateH", "IsingLogLikelihood", "IsingLogLikelihoodhJ", "IsingLogLikelihoodProbs", "binaryLogLikelihood"))
ejaList <- sapply(1:10, list)
res <- parSapply(cl, ejaList, bootstrapLikelihoodRatio, p = p, d = d, N = N, largeModel = largeModel, nullModel = nullModel, epsilon = epsilon, GLarge = GLarge, GNull = GNull)
stopCluster(cl)

bootstrapLikelihoodRatio(1, p, d, N, largeModel, nullModel, epsilon, GLarge, GNull)

