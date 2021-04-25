#Step 1. In case you don't have it yet, download R from http://www.r-project.org/. An excellent online introduction to R by Rob Cabacoff is available from http://www.statmethods.net/index.html. 

#Step 2. Install and load required packages and functions for correlations and partial correlations graphs by running the code below.
install.packages("qgraph") #see http://www.jstatsoft.org/v48/i04/ for details 
library(qgraph) #load qgraph

#Step 3. Read in NCS-R data (http://www.hcp.med.harvard.edu/ncs/). In these data, that are missing by design due to skip questions have been recoded to zero (indicating absence of symptom) and duration criteria have been removed.
setwd("D:/universitet/matematik/Speciale/R/data")
ncsdata=read.table(file="DepressionAnxiety.txt") #read in data
colnames(ncsdata)=c("depr", "inte", "weig", "mSle", "moto", "mFat", "repr", "conc", "suic", "anxi", "even", "ctrl", "edge", "gFat", "irri", "gCon", "musc", "gSle") #Define variable names.

#Step 4. Produce Figure 2a
n=15 #number of symptoms: 6 MD, 3 bridge symptoms and 6 GAD
labels=c("depr","inte","weig","moto","repr","suic","slee","conc","fati","anxi","even","ctrl","irri","musc","edge")
groups=list(MD=1:6,Bridge=7:9,GAD=10:15)
dataAdjacency=matrix(c( #adjacency matrix
0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,
1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,
1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,
1,1,1,0,1,1,1,1,1,0,0,0,0,0,0,
1,1,1,1,0,1,1,1,1,0,0,0,0,0,0,
1,1,1,1,1,0,1,1,1,0,0,0,0,0,0,
1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,
0,0,0,0,0,0,1,1,1,0,1,1,1,1,1,
0,0,0,0,0,0,1,1,1,1,0,1,1,1,1,
0,0,0,0,0,0,1,1,1,1,1,0,1,1,1,
0,0,0,0,0,0,1,1,1,1,1,1,0,1,1,
0,0,0,0,0,0,1,1,1,1,1,1,1,0,1,
0,0,0,0,0,0,1,1,1,1,1,1,1,1,0
),n,n)
qgraph(dataAdjacency,filename="Figure 2A",filetype="pdf",groups=groups,legend=T,labels=labels,color=c("red","green","lightblue"))

#Step 5. Produce the correlation network of Figure 2b.
ncscorrelationgraph=qgraph(cor(ncsdata),layout="spring",labels=colnames(ncsdata))

# own code:
selectVar <- !(attributes(ncsdata)$names == "anxi" | attributes(ncsdata)$names == "even")
ncsDataReduced <- subset(ncsdata, select = selectVar)
for(i in 1:16){
  for(j in 1:9282){
    if(ncsDataReduced[j, i] == 0) ncsDataReduced [j, i] <- -1
  }
}
ncsAvg <- apply(ncsDataReduced, 2, mean)
ncsM <- matrix(0, 16, 16)
for(i in seq_len(dim(ncsDataReduced)[1])) ncsM <- ncsM + as.numeric(ncsDataReduced[i, ]) %*% t(as.numeric(ncsDataReduced[i, ])) / dim(ncsDataReduced)[1]
ncsG <- list(1:16, 
             cbind(
               c(rep(1,15),rep(2,14),rep(3,13),rep(4,12),rep(5,11),rep(6,10),rep(7,9),rep(8,8),rep(9,7),rep(10,6),rep(11,5),rep(12,4),rep(13,3),rep(14,2),15),
               c(2:16,3:16,4:16,5:16,6:16,7:16,8:16,9:16,10:16,11:16,12:16,13:16,14:16,15:16,16)
             ))

#benchmark rcpp implementation:
microbenchmark(Ising_MLE_rcpp(xBar = ncsAvg, M = ncsM, G = ncsG, epsilon = 1e-9, epsilon2 = 1e-9, maxiter = 200),times = 10) # using strick convergence criteria
microbenchmark(Ising_MLE_rcpp(ncsAvg, ncsM, ncsG, 1e-4),times = 10)
microbenchmark(Ising_MLE_rcpp(xBar = ncsAvg, M = ncsM, G = ncsG, epsilon = 1e-4, epsilon2 = 1e-4),times = 100) # using critera reported in article
microbenchmark(Ising_MLE_rcpp(xBar = ncsAvg, M = ncsM, G = ncsG, epsilon = 7e-4, epsilon2 = 7e-4),times = 100) # converges in 28 iterations

# produce plot from figure 2 in article (using strict criteria):
IsingMLE <- Ising_MLE_rcpp(data = ncsDataReduced, G = ncsG, epsilon = 1e-9, epsilon2 = 1e-9, maxiter = 200)
J_hat <- IsingMLE$J
#colnames(J_hat) <- c("depr", "inte", "weig", "mSle", "moto", "mFat", "repr", "conc", "suic", "ctrl", "edge", "gFat", "irri", "gCon", "musc", "gSle")
#rownames(J_hat) <- c("depr", "inte", "weig", "mSle", "moto", "mFat", "repr", "conc", "suic", "ctrl", "edge", "gFat", "irri", "gCon", "musc", "gSle")
ncscorrelationgraph=qgraph(J_hat,layout="spring",labels=c("depr", "inte", "weig", "mSle", "moto", "mFat", "repr", "conc", "suic", "ctrl", "edge", "gFat", "irri", "gCon", "musc", "gSle"))
# produce plot from figure 2 in article (using criteria matching article):
IsingMLE <- Ising_MLE_rcpp(data = ncsDataReduced, G = ncsG, epsilon = 1e-4, epsilon2 = 1e-4)
J_hat <- IsingMLE$J
ncscorrelationgraph=qgraph(J_hat,layout="spring",labels=c("depr", "inte", "weig", "mSle", "moto", "mFat", "repr", "mCon", "suic", "ctrl", "edge", "gFat", "irri", "gCon", "musc", "gSle"))
# produce correlation network from figure 2 in article:
ncscorrelationgraph=qgraph(cor(ncsDataReduced),layout="spring",labels=c("depr", "inte", "weig", "mSle", "moto", "mFat", "repr", "mCon", "suic", "ctrl", "edge", "gFat", "irri", "gCon", "musc", "gSle"))

# not assuming MTP_2:
IsingMLEnoMTP2 <- Ising_MLE_noMTP2_rcpp(data = ncsDataReduced, G = ncsG, epsilon = 1e-9, epsilon2 = 1e-9, maxiter = 2000)
microbenchmark(Ising_MLE_noMTP2_rcpp(data = ncsDataReduced, G = ncsG, epsilon = 1e-9, epsilon2 = 1e-9, maxiter = 2000), times = 1)
microbenchmark(Ising_MLE_noMTP2_rcpp(data = ncsDataReduced, G = ncsG, epsilon = 1e-4, epsilon2 = 1e-4, maxiter = 2000), times = 1)
