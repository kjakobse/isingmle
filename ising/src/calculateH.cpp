#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the natural parameter h in the Ising model from the distribution p
//'
//' \code{calculateH} calculates the natural parameter h in the Ising model from the distribution p.
//'
//' @param d Integer specifying the number of binary variables in the sample.
//' @param p NumericVector specifying the distribution of the d binary variables.
//' @return \code{calculateH} returns a NumericVector with the natural parameter h.
//' @export
//'
// [[Rcpp::export]]
NumericVector calculateH(int d, NumericVector p) {
  NumericVector h(d);

  unsigned long long int p2powd = 1 << d;
  // loop over d and determine h for each i using that in the Ising model, it holds for any x with x_i=1 and y equal x up to the i'th coordinate, that log(p(x)p(-y)/(p(-x)p(y)))=4h_i.
  // x is chosen such that all coordinates other than x_i are -1.
  for (int i = 0; i < d; i++) {
    unsigned long long int p2powDMinusI = 1 << (d-i-1);
    h[i] = log(p[p2powDMinusI] * p[p2powd - 1] / (p[p2powd - p2powDMinusI - 1] * p[0])) / 4;
  }

  return(h);
}
