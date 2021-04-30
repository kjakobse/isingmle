#include <Rcpp.h>
using namespace Rcpp;

//' calculatexBar <to be documented>
//'
//' @param data <to be documented>
//' @param length <to be documented>
//' @param d <to be documented>
//' @return NumericVector <to be documented>
//' @export
//'
// [[Rcpp::export]]
NumericVector calculatexBar(NumericMatrix data, unsigned long long int length, int d) {
  NumericVector xBar(d);

  for (int t = 0; t < d; t++) {
    double dEntry = 0;

    for (unsigned long long int u = 0; u < length; u++) {
      dEntry += data(u, t);
    }
    xBar[t] = dEntry / length;
  }

  return (xBar);
}
