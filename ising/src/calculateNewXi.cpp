#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the second moment from the distribution
//'
//' \code{calculateNewXi} takes the distribution of d binary variables and
//' returns the matrix of second moments.
//'
//' @param p Numeric vector specifying the distribution of d binary variables.
//' @param d Integer specifying the number of binary variables.
//' @return \code{calculateNewXi} returns a numeric d by d matrix with the
//' second moments.
//' @export
//'
// [[Rcpp::export]]
NumericMatrix calculateNewXi(NumericVector p, int d) {
  NumericMatrix xi(d, d);
  unsigned long long int pLength = p.size();
  // loop over upper triangular part of Xi and calculate second moment:
  for (int t = 0; t < d; t++) {
    unsigned long long int tMask = 1 << (d - t - 1);
    for (int u = (t+1); u < d; u++) {
      unsigned long long int uMask = 1 << (d - u - 1);

      double entryValue = 0;
      // loop over entries in p and add if t'th and u'th observations are the
      // same, subtract if they are opposite signs:
      for (unsigned long long int v = 0; v < pLength; v++) {
        if (((bool)(v & tMask)) == ((bool)(v & uMask))) {
          entryValue += p[v]; // observations have same value
        } else {
          entryValue -= p[v]; // observations have opposite value
        }
      }
      xi(t, u) = entryValue;
      // save symmetric value:
      xi(u, t) = entryValue;

    }
  }
  // Diagonal elements are always 1:
  for (int s = 0; s < d; s++) {
    xi(s, s) = 1;
  }

  return(xi);
}
