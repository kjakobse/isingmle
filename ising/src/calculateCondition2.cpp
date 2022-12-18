#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the differences between M and Xi
//'
//' \code{calculateCondition2} calculates the differences between entries in
//' the matrices M and Xi given by E. This function is intended to be used by
//' IsingMLE and IsingMLEmtp2.
//'
//' @param e NumericMatrix containing in the rows the indices for the entries
//' in m and xi.
//' @param m NumericMatrix containing a sample second moment.
//' @param xi NumericMatrix to be compared with m.
//' @return \code{calculateCondition2} returns a vector with the difference
//' between entries in m and xi for those entries specified by e.
//' @export
//' @importFrom Rcpp evalCpp
//' @useDynLib ising, .registration = TRUE
//'
// [[Rcpp::export]]
NumericVector calculateCondition2(NumericMatrix e,
                                  NumericMatrix m,
                                  NumericMatrix xi) {
  int nRows = e.nrow();

  NumericVector condition2(nRows);
// For each row in e, calculate the difference between the given entries of m
// and xi. Note that the entries in e are assumed 1-indexed, so 1 is subtracted
// as c++ uses 0-indexing:
  for (int i = 0; i < nRows; i++) {
    condition2[i] = m(e(i, 0) - 1, e(i, 1) - 1) - xi(e(i, 0) - 1, e(i, 1) - 1);
  }

  return(condition2);
}
