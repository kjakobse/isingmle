#include <Rcpp.h>
using namespace Rcpp;

//' Calculate empirical marginal distributions of pairs in ePlus
//'
//' \code{CalculateEmpirical} calculates the empirical distribution of pairs of binary variables specified by ePlus from the sample first and second moments.
//'
//' ePlus should contain the variable pairs for which the empirical distribution is to be calculated. They should be coded as integers from 1 up to the number of binary variables.
//' The sample moments should be ordered according to this coding, so the sample mean of variable 1 should have index 1 in xBar (using 1-indexing).
//'
//' @param ePlus NumericMatrix containing pairs of variables in the rows.
//' @param m NumericMatrix containing the sample second moment.
//' @param xBar NumericVector containing the sample first moment.
//' @return \code{CalculateEmpirical} returns a List containing a NumericMatrix with the calculated empirical distributions and an integer which is 1 of any of the empirical distributions contain zeroes and is 0 otherwise.
//' @export
//'
// [[Rcpp::export]]
List calculateEmpirical(NumericMatrix ePlus, NumericMatrix m, NumericVector xBar) {
// determine the number of pairs to calculate the empirical distribution of:
  int ePlusDim = ePlus.nrow();
// initiate the integer indicating if any empirical distributions contain zeroes:
  int containsZeroes = 0;

  NumericMatrix empirical(ePlusDim, 4);
  for (int t = 0; t < ePlusDim; t++) {
    // for each row in ePlus, extract the variable pair, and calculate the empirical distribution from the corresponding entries in the sample moments:
    int i = ePlus(t, 0) - 1;
    int j = ePlus(t, 1) - 1;
    empirical(t, 0) = (1 + xBar[i] + xBar[j] + m(i, j)) / 4;
    empirical(t, 1) = (1 + xBar[i] - xBar[j] - m(i, j)) / 4;
    empirical(t, 2) = (1 - xBar[i] + xBar[j] - m(i, j)) / 4;
    empirical(t, 3) = (1 - xBar[i] - xBar[j] + m(i, j)) / 4;
    // check if the empirical distribution of pair ij contains zeroes and set containsZeroes to one if it does:
    if (abs(empirical(t, 0)) <= 1e-15 || abs(empirical(t, 1)) <= 1e-15 || abs(empirical(t, 2)) <= 1e-15 || abs(empirical(t, 3)) <= 1e-15) {
      containsZeroes = 1;
    }
    if (empirical(t, 0) < -1e-15 || empirical(t, 1) < -1e-15 || empirical(t, 2) < -1e-15 || empirical(t, 3) < -1e-15) {
      NumericMatrix error = NumericMatrix(0, 0);
      return List::create(_["empirical"] = error);
    }
  }
  return List::create(_["empirical"] = empirical, _["containsZeroes"] = containsZeroes);
}
