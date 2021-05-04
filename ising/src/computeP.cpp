#include <Rcpp.h>
using namespace Rcpp;

//' Calculate distribution from canonical parameters
//'
//' \code{computeP} calculates the distribution of an Ising model from the canonical parameters h and J.
//'
//' The distribution p is returned in a vector where the index (using 0-indexing) corresponds to the base 10 value of the binary observation. \cr
//' The probability p(x) is obtained using the formula \eqn{p(x) = e^(h^T * x + x^T * J * x / 2 - A(h,J))}, where A(h,J) is a normalisation constant.
//'
//' @param h Numeric vector specifying the canonical parameter h.
//' @param J Matrix specifying the canonical parameter J.
//' @return \code{computeP} returns a numeric vector with the distribution of the Ising model with canonical parameters h and J.
//' @export
//'
// [[Rcpp::export]]
NumericVector computeP(NumericVector h, NumericMatrix J) {
  // initialise number of variables and number of possible observations from size of h:
  int d = h.size();
  unsigned long long int pLength = pow(2, d);
  NumericVector p(pLength);
  double pSum = 0;
  // loop over all possible observations:
  for (unsigned long long int t = 0; t < pLength; t++) {
    double pEntry = 0;
    // loop over indices in h and first indices in J:
    for (int u = 0; u < d; u++) {
      // mask for the index u. << is the bitwise and operator and shifts the bits to the left:
      unsigned long long int tMaskU = 1 << u;
      // Determine the value of the u'th binary variable for the current index of p:
      int observationValueU = -1;
      if (t & tMaskU) {
        observationValueU = 1;
      }
      // loop over second indices in J:
      for(int v = 0; v < d; v++) {
        // mask for index v then determine the value of the v'th binary variable for the current index of p:
        unsigned long long int tMaskV = 1 << v;
        int observationValueV = -1;
        if (t & tMaskV) {
          observationValueV = 1;
        }
        // add the J_uv term to p:
        pEntry += 0.5 * J(d-u-1, d-v-1) * observationValueU * observationValueV;
      }
      // add the h_u term to p
      pEntry += h(d-u-1) * observationValueU;
    }
    // take the exponential of pEntry:
    p[t] = exp(pEntry);
    // add the term to the normalisation constant:
    pSum += exp(pEntry);
  }
  // normalise the entries in p to sum to 1:
  for (unsigned long long int t = 0; t < pLength; t++) {
    p[t] = p[t] / pSum;
  }

  return(p);
}
