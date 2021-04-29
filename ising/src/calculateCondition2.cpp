#include <Rcpp.h>
using namespace Rcpp;

//' calculateCondition2 <to be documented>
//'
//' @param e <to be documented>
//' @param m <to be documented>
//' @param xi <to be documented>
//' @return NumericVector <to be documented>
//' @export
//
// [[Rcpp::export]]
NumericVector calculateCondition2(NumericMatrix e, NumericMatrix m, NumericMatrix xi) {
  int nRows = e.nrow();

  NumericVector condition2(nRows);
  for (int i = 0; i < nRows; i++) {
    condition2[i] = m(e(i, 0) - 1, e(i, 1) - 1) - xi(e(i, 0) - 1, e(i, 1) - 1);
  }

  return(condition2);
}
