#include <Rcpp.h>
using namespace Rcpp;

//' calculateNewXi <to be documented>
//'
//' @param p <to be documented>
//' @param d <to be documented>
//' @return NumericMatrix <to be documented>
//' @export
//
// [[Rcpp::export]]
NumericMatrix calculateNewXi(NumericVector p, int d) {
  NumericMatrix xi(d, d);
  unsigned long long int pLength = p.size();

  for (int t = 0; t < d; t++) {
    unsigned long long int tMask = 1 << (d - t - 1);
    for (int u = t; u < d; u++) {
      unsigned long long int uMask = 1 << (d - u - 1);

      double entryValue = 0;
      for (unsigned long long int v = 0; v < pLength; v++) {
        if (((bool)(v & tMask)) == ((bool)(v & uMask))) {
          entryValue += p[v]; // observations have same value
        } else {
          entryValue -= p[v]; // observations have opposite value
        }
      }
      xi(t, u) = entryValue;
      if (t != u) { // save symmetric value
        xi(u, t) = entryValue;
      }
    }
  }

  return(xi);
}
