#include <Rcpp.h>
using namespace Rcpp;

//' calculateCondition <to be documented>
//'
//' @param eHatOmitNA <to be documented>
//' @param m <to be documented>
//' @param xi <to be documented>
//' @return NumericVector <to be documented>
//' @export
//
// [[Rcpp::export]]
NumericVector calculateCondition(NumericMatrix eHatOmitNA, NumericMatrix m, NumericMatrix xi) {
  int nRows = eHatOmitNA.nrow();

  NumericVector condition(nRows);
  for (int i = 0; i < nRows; i++) {
    condition[i] = abs(m(eHatOmitNA(i, 0) - 1, eHatOmitNA(i, 1) - 1) - xi(eHatOmitNA(i, 0) - 1, eHatOmitNA(i, 1) - 1));
  }

  return(condition);
}
