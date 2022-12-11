#include <Rcpp.h>
using namespace Rcpp;

//' Create initial distribution for MLE algorithm
//'
//' \code{createP} is intended to be used with \code{isingMLE} and \code{isingMLEmtp2}. It creates the vector p with the initial distribution that the algorithms start from.
//'
//' The initial distribution has the form \eqn{p(x) = 2^{-|V|}\prod_{v \in V}(1 + x_v \mu_v)}.
//'
//' @param mu Numeric vector specifying the mean of the distribution.
//' @param d Integer specifying the number of binary variables.
//' @return \code{createP} returns a numeric vector with the distribution p.
//' @export
//'
// [[Rcpp::export]]
NumericVector createP(NumericVector mu, int d) {
  unsigned long long int pLength = pow(2, d);
  NumericVector p(pLength);
  // Initialise the normalisation constant:
  double normalisationConstant = pow(2, -d);
  // calculate the probability for each possible outcome:
  for (unsigned long long int t = 0; t < pLength; t++) {
    double pEntry = 1;
    // loop over the binary variables and determine their value for the current observation:
    for (int u = 0; u < d; u++) {
      // mask for the index u. << is the bitwise and operator and shifts the bits to the left:
      unsigned long long int tMask = 1 << u;
      // use the mask to determine the value of the u'th binary variable:
      int observationValue = -1;
      if (t & tMask) {
        observationValue = 1;
      }
      // add term for the u'th binary variable to the product:
      double value = 1 + observationValue * mu[d - u - 1];
      pEntry *= value;
    }
    // normalise the distribution:
    p[t] = normalisationConstant * pEntry;
  }

  return(p);
}
