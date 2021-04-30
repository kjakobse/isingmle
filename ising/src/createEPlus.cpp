#include <Rcpp.h>
using namespace Rcpp;

//' createEPlus <to be documented>
//'
//' @param e <to be documented>
//' @param m <to be documented>
//' @param xBar <to be documented>
//' @return NumericMatrix <to be documented>
//' @export
//'
// [[Rcpp::export]]
NumericMatrix createEPlus(NumericMatrix e, NumericMatrix m, NumericVector xBar) {
  int nRows = e.nrow();
  int nCols = e.ncol();

  NumericVector includeRows(nRows);
  int nIncludedRows = 0;
  for (int i = 0; i < nRows; i++) {
    if (m(e(i, 0) - 1, e(i, 1) - 1) > (xBar[e(i, 0) - 1] * xBar[e(i, 1) - 1])) {
      includeRows[i] = 1;
      nIncludedRows++;
    }
  }
  NumericMatrix ePlus(nIncludedRows, nCols);
  int k = 0;
  for (int i = 0; i < nRows; i++) {
    if (includeRows[i] == 1) {
      for (int j = 0; j < nCols; j++) {
        ePlus(k, j) = e(i, j);
      }
      k++;
    }
  }

  return(ePlus);
}
