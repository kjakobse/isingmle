library(Rcpp)

# c++ functions for the IsingMLE function:

cppFunction('
  List calculateNewPAndEHatnoMTP2(NumericMatrix ePlus, int d, NumericMatrix e, NumericVector p, NumericMatrix eHat, NumericMatrix capitalJ, int iter) {
    int ePlusDim = ePlus.nrow();
    unsigned long long int pLength = p.size();
    for (int t = 0; t < ePlusDim; t++) {
      int i = ePlus(t, 0);
      int j = ePlus(t, 1);
      int iCInternal = d - i;
      int jCInternal = d - j;
      unsigned long long int iMask = 1 << iCInternal; 
      unsigned long long int jMask = 1 << jCInternal; 

      double p00 = 0.0;
      double p01 = 0.0;
      double p10 = 0.0;
      double p11 = 0.0;
      for (unsigned long long int u = 0; u < pLength; u++) {
        if (u & iMask) {
          if (u & jMask) {
            p11 += p[u];
          } else {
            p10 += p[u];
          }
        } else {
          if (u & jMask) {
            p01 += p[u];
          } else {
            p00 += p[u];
          }
        }
      }

      double q11 = e(t, 0) / p11;
      double q10 = e(t, 1) / p10;
      double q01 = e(t, 2) / p01;
      double q00 = e(t, 3) / p00;
      
      unsigned long long int v1 = 0; 
      unsigned long long int v2 = iMask; 
      unsigned long long int v3 = jMask; 
      unsigned long long int v4 = iMask | jMask; 
      if(iter == 0) {
        capitalJ(i-1, j-1) = 0.25 * log(p[v1] * p[v4] / (p[v2] * p[v3]));
        capitalJ(j-1, i-1) = 0.25 * log(p[v1] * p[v4] / (p[v2] * p[v3]));
      }

      for (unsigned long long int u = 0; u < pLength; u++) {
        if (u & iMask) {
          if (u & jMask) {
            p[u] = p[u] * q11;
          } else {
            p[u] = p[u] * q10;
          }
        } else {
          if (u & jMask) {
            p[u] = p[u] * q01;
          } else {
            p[u] = p[u] * q00;
          }
        }
      }
      
      capitalJ(i-1, j-1) = 0.25 * log(p[v1] * p[v4] / (p[v2] * p[v3]));
      capitalJ(j-1, i-1) = 0.25 * log(p[v1] * p[v4] / (p[v2] * p[v3]));
      
      if(abs(capitalJ(i-1, j-1)) < 1e-15) {
        eHat(t, 0) = NA_REAL;
        eHat(t, 1) = NA_REAL;
      } else {
        eHat(t, 0) = i;
        eHat(t, 1) = j;
      }
    }
    
    return List::create(_["p"] = p, _["eHat"] = eHat, _["J"] = capitalJ);
  }
')

cppFunction('
  List calculateNewPAndEHatnoMTP2Boundary(NumericMatrix ePlus, int d, NumericMatrix e, NumericVector p, NumericMatrix eHat) {
    int ePlusDim = ePlus.nrow();
    unsigned long long int pLength = p.size();
    for (int t = 0; t < ePlusDim; t++) {
      int i = ePlus(t, 0);
      int j = ePlus(t, 1);
      int iCInternal = d - i;
      int jCInternal = d - j;
      unsigned long long int iMask = 1 << iCInternal; 
      unsigned long long int jMask = 1 << jCInternal; 

      double p00 = 0.0;
      double p01 = 0.0;
      double p10 = 0.0;
      double p11 = 0.0;
      for (unsigned long long int u = 0; u < pLength; u++) {
        if (u & iMask) {
          if (u & jMask) {
            p11 += p[u];
          } else {
            p10 += p[u];
          }
        } else {
          if (u & jMask) {
            p01 += p[u];
          } else {
            p00 += p[u];
          }
        }
      }
      
      double q11;
      double q10;
      double q01;
      double q00;
      
      if(abs(e(t, 0) < 1e-15) || abs(p11) < 1e-15) {
        q11 = 0;
      } else {
        q11 = e(t, 0) / p11;
      }
      if(abs(e(t, 1) < 1e-15) || abs(p10) < 1e-15) {
        q10 = 0;
      } else {
        q10 = e(t, 1) / p10;
      }
      if(abs(e(t, 2) < 1e-15) || abs(p01) < 1e-15) {
        q01 = 0;
      } else {
        q01 = e(t, 2) / p01;
      }
      if(abs(e(t, 3) < 1e-15) || abs(p00) < 1e-15) {
        q00 = 0;
      } else {
        q00 = e(t, 3) / p00;
      }
      eHat(t, 0) = i;
      eHat(t, 1) = j;

      for (unsigned long long int u = 0; u < pLength; u++) {
        if (u & iMask) {
          if (u & jMask) {
            p[u] = p[u] * q11;
          } else {
            p[u] = p[u] * q10;
          }
        } else {
          if (u & jMask) {
            p[u] = p[u] * q01;
          } else {
            p[u] = p[u] * q00;
          }
        }
      }
    }
    
    return List::create(_["p"] = p, _["eHat"] = eHat);
  }
')

# Function fitting without MTP2 condition:
IsingMLE <- function(G, xBar = NULL, M = NULL, data = NULL, epsilon = 1e-4, epsilon2 = 1e-4, maxiter = 100){
  V <- G[[1]]
  E <- G[[2]]
  
  if(is.null(data)) {
    if(is.null(xBar) | is.null(M)) {
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
  for(i in (1+seq_len(d-1))) {
    for(j in seq_len(i-1))
      mCond[i] <- M[i, j]
  }
  #if(!(all(abs(mCond) < 1))) {
  #  stop(cat("input doesn't fulfill conditions for existance of MLE.\nEmpirical correlation matrix contains 1 or -1."))
  #}
  
  Xi <- diag(d)
  
  p <- createP(mu, d)
  
  eHat <- matrix(NA, nrow(E), 2)
  eHatOmitNA <- na.omit(eHat)
  
  condition <- c(Inf)
  
  empiricalList <- CalculateEmpirical(E, M, xBar)
  empirical <- empiricalList$empirical
  containsZeroes <- empiricalList$containsZeroes
  if (ncol(empirical) == 0) {
    stop(cat("input doesn't fulfill conditions for existance of MLE.\n Produced negative values of the empirical distribution"))
  }
  
  if(containsZeroes) {
    iter <- 0L
    while((iter < maxiter) & (max(abs(mu-xBar)) > epsilon | max(condition) > epsilon)){ 
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
    
    if(iter >= maxiter) warning("MLE did not converge; maximum number of iterations reached")
    warning("Maximum likelihood estimate is on the boundary")
    
    return(list(p_hat = p, G_hat = list(V_hat = V, E_hat = eHatOmitNA), mu_hat = mu, Xi_hat = Xi, number_of_iterations = iter))
  } else {
    iter <- 0L
    capitalJ <- matrix(0, d, d)
    while((iter < maxiter) & (max(abs(mu-xBar)) > epsilon | max(condition) > epsilon)){ 
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
  
    if(iter >= maxiter) warning("MLE did not converge; maximum number of iterations reached")
  
    return(list(p_hat = p, G_hat = list(V_hat = V, E_hat = eHatOmitNA), mu_hat = mu, Xi_hat = Xi, h_hat = h, J_hat = capitalJ, number_of_iterations = iter))
  }
}


