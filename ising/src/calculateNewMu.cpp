#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the mean vector from the distribution
//'
//' \code{calculateNewMu} calculates the mean for each binary variable from the
//' distribution of d binary variables.
//'
//' @param p Numeric vector specifying the distribution of d binary variables.
//' @param d Integer specifying the number of binary variables.
//' @return \code{calculateNewMu} returns a numeric vector with the means.
//' @author Kim Daniel Jakobsen
//' @examples
//' 1+1
//'
//' @export
//'
// [[Rcpp::export]]

NumericVector calculateNewMu(NumericVector p, int d) {
  NumericVector mu(d);
  unsigned long long int pLength = p.size();
// run through each variable and add up observation*probability:
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
