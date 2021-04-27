library(Rcpp)

# c++ function for calculating the likelihood: IsingLogLikelihood(matrix(0,1,1), h, J)

cppFunction('
  double IsingLogLikelihoodhJ(NumericMatrix data, NumericVector h, NumericMatrix J) {
    int d = h.size();
    unsigned long long int nrow = data.nrow();
    unsigned long long int pLength = pow(2, d);
    double pSum = 0;
    double logLikelihood = 0;

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
      pSum += exp(pEntry);
    }

    double A = log(pSum);

    for(unsigned long long int t = 0; t < nrow; t++) {
      double dataEntry = 0;
      for (int u = 0; u < d; u++) {

        for(int v = 0; v < d; v++) {

          dataEntry += 0.5 * J(u, v) * data(t, u) * data(t, v);
        }
        dataEntry += h(u) * data(t, u);
      }
      dataEntry -= A;
      logLikelihood += dataEntry;
    }

    return(logLikelihood);
  }
')

cppFunction('
  double IsingLogLikelihoodProbs(NumericMatrix data, NumericVector p) {
    int d = data.ncol();
    unsigned long long int nrow = data.nrow();
    double logLikelihood = 0;

    for(unsigned long long int t = 0; t < nrow; t++) {
      double pEntry = 0;
      for (int u = 0; u < d; u++) {
          if(data(t, u) == 1) {
            pEntry += pow(2, d-u-1);
        }
      }
      logLikelihood += log(p[pEntry]);
    }

    return(logLikelihood);
  }
')

IsingLogLikelihood <- function(data, h = NULL, J = NULL, p = NULL) {
  if(is.null(h) & is.null(J) & !is.null(p)) {
    IsingLogLikelihood <- IsingLogLikelihoodProbs(data, p)
  } else if(!is.null(h) & !is.null(J) & is.null(p)) {
    IsingLogLikelihood <- IsingLogLikelihoodhJ(data, h, J)
  } else {
    stop("Must specify either the natural parameters h and J or the probability vector p")
  }
  IsingLogLikelihood
}
