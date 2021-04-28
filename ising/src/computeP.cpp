#include <Rcpp.h>
using namespace Rcpp;

//' computeP <to be documented>
//'
//' @param h <to be documented>
//' @param J <to be documented>
//' @return <to be documented>
//' @export
//
// [[Rcpp::export]]
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
