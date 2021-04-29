#include <Rcpp.h>
using namespace Rcpp;

//' isingLogLikelihoodProbs <to be documented>
//'
//' @param data <to be documented>
//' @param p <to be documented>
//' @return double <to be documented>
//' @export
//
// [[Rcpp::export]]
double isingLogLikelihoodProbs(NumericMatrix data, NumericVector p) {
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
