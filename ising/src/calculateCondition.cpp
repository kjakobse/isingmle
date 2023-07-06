#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the distances between M and Xi
//'
//' \code{calculateCondition} calculates the absolute differences between
//' entries in the matrices M and Xi given by eHatOmitNA. This function is
//' intended to be used by IsingMLE and IsingMLEmtp2.
//'
//' @param eHatOmitNA NumericMatrix containing in the rows the indices for the
//' entries in m and xi.
//' @param m NumericMatrix containing a sample second moment.
//' @param xi NumericMatrix to be compared with m.
//' @return \code{calculateCondition} returns a vector with the distance
//' between entries in m and xi for those entries specified by eHatOmitNA.
//' @export
//' @author Kim Daniel Jakobsen
//' @examples
//' 1+1
//'
// [[Rcpp::export]]

NumericVector calculateCondition(NumericMatrix eHatOmitNA,
                                 NumericMatrix m,
                                 NumericMatrix xi) {
  int nRows = eHatOmitNA.nrow();

  NumericVector condition(nRows);
  // For each row in eHatOmitNA, calculate the distance between the given
  // entries of m and xi. Note that the entries in eHatOmitNA are assumed
  // 1-indexed, so 1 is subtracted as c++ uses 0-indexing:
  for (int i = 0; i < nRows; i++) {
    condition[i] = abs(
      m(eHatOmitNA(i, 0) - 1, eHatOmitNA(i, 1) - 1) -
        xi(eHatOmitNA(i, 0) - 1, eHatOmitNA(i, 1) - 1)
    );
  }

  return(condition);
}
