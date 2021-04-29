#include <Rcpp.h>
using namespace Rcpp;

//' createP <to be documented>
//'
//' @param mu <to be documented>
//' @param d <to be documented>
//' @return NumericVector <to be documented>
//' @export
//
// [[Rcpp::export]]
NumericVector createP(NumericVector mu, int d) {
  unsigned long long int pLength = pow(2, d);
  NumericVector p(pLength);
  double normalisationConstant = pow(2, -d);
  for (unsigned long long int t = 0; t < pLength; t++) {
    double pEntry = 1;
    for (int u = 0; u < d; u++) {
      unsigned long long int tMask = 1 << u;
      int observationValue = -1;
      if (t & tMask) {
        observationValue = 1;
      }
      double value = 1 + observationValue * mu[d - u - 1];
      pEntry *= value;
    }
    p[t] = normalisationConstant * pEntry;
  }

  return(p);
}
