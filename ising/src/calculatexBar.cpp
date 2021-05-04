#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the sample first moment
//'
//' \code{calculatexBar} calculates the sample first moment for a data set with samples from d binary variables.
//'
//' @param data Matrix containing samples from d binary variables in the rows.
//' @param length Integer specifying the number of samples in the data set.
//' @param d Integer specifying the number of binary variables.
//' @return \code{calculatexBar} returns a numeric vector with the sample first moment.
//' @export
//'
// [[Rcpp::export]]
NumericVector calculatexBar(NumericMatrix data, unsigned long long int length, int d) {
  NumericVector xBar(d);
  // For each binary variable, run through the samples, add them up, and divide by the number of samples:
  for (int t = 0; t < d; t++) {
    double dEntry = 0;

    for (unsigned long long int u = 0; u < length; u++) {
      dEntry += data(u, t);
    }
    xBar[t] = dEntry / length;
  }

  return (xBar);
}
