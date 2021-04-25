library(Rcpp)

IntToSign <- function(x, dig) {
  i <- 0L
  string <- rep(-1, dig)
  while (x > 0) {
    if(x %% 2L == 0) string[dig - i] <- -1
    else string[dig - i] <- 1
    x <- x %/% 2L
    i <- i + 1L
  }
  string
}

cppFunction('
  NumericVector computeP(NumericVector h, NumericMatrix J) {
    int d = h.size();
    unsigned long long int pLength = pow(2, d);
    NumericVector p(pLength);
    double pSum = 0;

    for (unsigned long long int t = 0; t < pLength; t++) {
      double pEntry = 0;
      for (int u = 0; u < d; u++) {
        unsigned long long int tMaskU = 1 << u;
        int observationValueU = -1;
        if (t & tMaskU) {
          observationValueU = 1;
        }
        for(int v = 0; v < d; v++) {
          unsigned long long int tMaskV = 1 << v;
          int observationValueV = -1;
          if (t & tMaskV) {
            observationValueV = 1;
          }
          pEntry += 0.5 * J(d-u-1, d-v-1) * observationValueU * observationValueV;
        }
        pEntry += h(d-u-1) * observationValueU;
      }
      p[t] = exp(pEntry);
      pSum += exp(pEntry);
    }
    
    for (unsigned long long int t = 0; t < pLength; t++) {
      p[t] = p[t] / pSum;
    }
    
    return(p);
  }
')

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
