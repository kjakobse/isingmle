#include <Rcpp.h>
using namespace Rcpp;

//' calculatexBar <to be documented>
//'
//' @param data <to be documented>
//' @param length <to be documented>
//' @param d <to be documented>
//' @return NumericMatrix <to be documented>
//' @export
//'
// [[Rcpp::export]]
NumericMatrix calculateM(NumericMatrix data, unsigned long long int length, int d) {
  NumericMatrix matrixM(d, d);

  for (int t = 0; t < d; t++) {
    for (int u = 0; u < d; u++) {
      double entryValue = 0;
      for (unsigned long long int v = 0; v < length; v++) {
        //entryValue += data(v, t) * data(v, u);
        if (data(v, t) == data(v, u)) {
          entryValue++;
        } else {
          entryValue--;
        }
      }
      matrixM(t, u) = entryValue / length;
    }
  }

  return(matrixM);
}
