#include <Rcpp.h>
using namespace Rcpp;

//' Creates the matrix EPlus
//'
//' \code{createEPlus} is intended to be used with \code{isingMLE} and
//' \code{isingMLEmtp2}. It creates the matrix EPlus of edges uv where
//' \eqn{M_uv > \bar{x}_u \bar{x}_v}.
//'
//' @param e Matrix containing the edges of the graph G.
//' @param m Matrix containing the sample second moment.
//' @param xBar Numeric vector containing the sample first moment.
//' @return \code{createEPlus} returns the matrix EPlus.
//' @export
//'
// [[Rcpp::export]]
NumericMatrix createEPlus(NumericMatrix e,
                          NumericMatrix m,
                          NumericVector xBar) {
  int nRows = e.nrow();
  int nCols = e.ncol();

  // Initialise indicator of which rows to include in EPlus and a counter of
  // how many total:
  NumericVector includeRows(nRows);
  int nIncludedRows = 0;
  // Check condition for being in EPlus. If true, set includeRows to 1 and
  // increment total row counter:
  for (int i = 0; i < nRows; i++) {
    if (m(e(i, 0) - 1, e(i, 1) - 1) > (xBar[e(i, 0) - 1] * xBar[e(i, 1) - 1])) {
      includeRows[i] = 1;
      nIncludedRows++;
    }
  }
  NumericMatrix ePlus(nIncludedRows, nCols);
  int k = 0;
  // Run through rows of e and add to EPlus if indicator includeRows is 1:
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
