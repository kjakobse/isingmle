#include <Rcpp.h>
using namespace Rcpp;

//' Calculate empirical marginal distributions of pairs in ePlus adding small epsilon to zeroes.
//'
//' \code{CalculateEmpiricalReplaceZeroes} calculates the empirical distribution of pairs of binary variables specified by ePlus from the sample first and second moments and adds a small epsilon if a probability is zero.
//'
//' ePlus should contain the variable pairs for which the empirical distribution is to be calculated. They should be coded as integers from 1 up to the number of binary variables.
//' The sample moments should be ordered according to this coding, so the sample mean of variable 1 should have index 1 in xBar (using 1-indexing).
//'
//' @param ePlus Matrix containing pairs of variables in the rows.
//' @param m Matrix containing the sample second moment.
//' @param xBar Numeric vector containing the sample first moment.
//' @param eps Double to replace zeroes with. Must be >0.
//' @return \code{CalculateEmpirical} returns a matrix with the calculated empirical distributions.
//' @export
//'
// [[Rcpp::export]]
NumericMatrix calculateEmpiricalReplaceZeroes(NumericMatrix ePlus, NumericMatrix m, NumericVector xBar, double eps) {
  int ePlusDim = ePlus.nrow();
  NumericMatrix empirical(ePlusDim, 4);
  for (int t = 0; t < ePlusDim; t++) {
    int i = ePlus(t, 0) - 1;
    int j = ePlus(t, 1) - 1;
    empirical(t, 0) = (1 + xBar[i] + xBar[j] + m(i, j)) / 4;
    empirical(t, 1) = (1 + xBar[i] - xBar[j] - m(i, j)) / 4;
    empirical(t, 2) = (1 - xBar[i] + xBar[j] - m(i, j)) / 4;
    empirical(t, 3) = (1 - xBar[i] - xBar[j] + m(i, j)) / 4;
    if (empirical(t, 0) <= 0) {
      empirical(t, 0) = eps;
    }
    if (empirical(t, 1) <= 0) {
      empirical(t, 1) = eps;
    }
    if (empirical(t, 2) <= 0) {
      empirical(t, 2) = eps;
    }
    if (empirical(t, 3) <= 0) {
      empirical(t, 3) = eps;
    }
  }
  return(empirical);
}
