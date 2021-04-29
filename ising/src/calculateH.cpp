#include <Rcpp.h>
using namespace Rcpp;

//' calculateH <to be documented>
//'
//' @param d <to be documented>
//' @param p <to be documented>
//' @return NumericVector <to be documented>
//' @export
//
// [[Rcpp::export]]
NumericVector calculateH(int d, NumericVector p) {
  NumericVector h(d);

  unsigned long long int p2powd = 1 << d;
  for (int i = 0; i < d; i++) {
    unsigned long long int p2powDMinusI = 1 << (d-i-1);
    h[i] = log(p[p2powDMinusI] * p[p2powd - 1] / (p[p2powd - p2powDMinusI - 1] * p[0])) / 4;
  }

  return(h);
}
