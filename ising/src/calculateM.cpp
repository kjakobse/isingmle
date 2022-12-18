#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the sample second moment
//'
//' \code{calculateM} calculates the sample second moment for a data set with
//' samples from d binary variables.
//'
//' @param data Matrix containing samples from d binary variables in the rows.
//' @param length Integer specifying the number of samples in the data set.
//' @param d Integer specifying the number of binary variables.
//' @return \code{calculateM} returns a matrix with the sample second moment.
//' @export
//'
// [[Rcpp::export]]
NumericMatrix calculateM(NumericMatrix data,
                         unsigned long long int length,
                         int d) {
  NumericMatrix matrixM(d, d);

// run through each entry in matrixM and for each check if each sample has t'th
// and u'th coordinate equal. If so add 1/length since the product of the two
// coordinates is 1.
// If not subtract 1/length since the product of the two coordinates is -1:
  for (int t = 0; t < d; t++) {
    for (int u = 0; u < d; u++) {
      double entryValue = 0;
      for (unsigned long long int v = 0; v < length; v++) {
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
