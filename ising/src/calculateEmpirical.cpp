#include <Rcpp.h>
using namespace Rcpp;

//' CalculateEmpirical <to be documented>
//'
//' @param ePlus <to be documented>
//' @param m <to be documented>
//' @param xBar <to be documented>
//' @return List <to be documented>
//' @export
//'
// [[Rcpp::export]]
List calculateEmpirical(NumericMatrix ePlus, NumericMatrix m, NumericVector xBar) {
  int ePlusDim = ePlus.nrow();
  int containsZeroes = 0;

  NumericMatrix empirical(ePlusDim, 4);
  for (int t = 0; t < ePlusDim; t++) {
    int i = ePlus(t, 0) - 1;
    int j = ePlus(t, 1) - 1;
    empirical(t, 0) = (1 + xBar[i] + xBar[j] + m(i, j)) / 4;
    empirical(t, 1) = (1 + xBar[i] - xBar[j] - m(i, j)) / 4;
    empirical(t, 2) = (1 - xBar[i] + xBar[j] - m(i, j)) / 4;
    empirical(t, 3) = (1 - xBar[i] - xBar[j] + m(i, j)) / 4;
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
