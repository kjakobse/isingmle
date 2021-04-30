#include <Rcpp.h>
using namespace Rcpp;

//' calculateNewMu <to be documented>
//'
//' @param p <to be documented>
//' @param d <to be documented>
//' @return NumericVector <to be documented>
//' @export
//'
// [[Rcpp::export]]
NumericVector calculateNewMu(NumericVector p, int d) {
  NumericVector mu(d);
  unsigned long long int pLength = p.size();

  for (int t = 0; t < d; t++) {
    unsigned long long int uMask = 1 << (d - t - 1);

    double dEntry = 0;
    for (unsigned long long int u = 0; u < pLength; u++) {
      if (u & uMask) {
        dEntry += p[u]; // Observation is +1 add p[u]
      } else {
        dEntry -= p[u]; // observation is -1 subtract p[u]
      }
    }
    mu[t] = dEntry;
  }

  return(mu);
}
