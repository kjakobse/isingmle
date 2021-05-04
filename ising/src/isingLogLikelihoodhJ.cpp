#include <Rcpp.h>
using namespace Rcpp;

//' Calculate loglikelihood for data set given canonical parameters
//'
//' \code{isingLogLikelihoodhJ} calculates the value of the loglikelihood function for a given data set in the natural parameters h and J.
//'
//' The loglikelihood is calculated using the formula \deqn{l(h, J\mid data) = \sum_{x}\exp(h^T x+x^T J x/2-A(h,J))}, where A(h,J) is a normalisation constant.
//'
//' @param data Matrix containing samples from d binary variables in the rows.
//' @param h Numeric vector specifying the canonical parameter h.
//' @param J Matrix specifying the canonical parameter J.
//' @return  \code{isingLogLikelihoodhJ} returns the numeric value of the loglikelihood function.
//' @export
//'
// [[Rcpp::export]]
double isingLogLikelihoodhJ(NumericMatrix data, NumericVector h, NumericMatrix J) {
  int d = h.size();
  unsigned long long int nrow = data.nrow();
  unsigned long long int pLength = pow(2, d);
  double pSum = 0;
  double logLikelihood = 0;
  // Calculate the normalisation constant A(h,J)by calculating the unnormalised distribution and summing:
  for (unsigned long long int t = 0; t < pLength; t++) {
    double pEntry = 0;
    for (int u = 0; u < d; u++) {
      // mask for the index u. << is the bitwise and operator and shifts the bits to the left:
      unsigned long long int tMaskU = 1 << u;
      // use the mask to determine the value of the u'th binary variable:
      int observationValueU = -1;
      if (t & tMaskU) {
        observationValueU = 1;
      }
      for(int v = 0; v < d; v++) {
        // mask for the index v.
        unsigned long long int tMaskV = 1 << v;
        // use the mask to determine the value of the v'th binary variable:
        int observationValueV = -1;
        if (t & tMaskV) {
          observationValueV = 1;
        }
        // add the J_uv term to the current entry in the unnormalised distribution:
        pEntry += 0.5 * J(d-u-1, d-v-1) * observationValueU * observationValueV;
      }
      // add the h_u term to the current entry in the unnormalised distribution:
      pEntry += h(d-u-1) * observationValueU;
    }
    // sum over the unnormalised distribution:
    pSum += exp(pEntry);
  }
  // A(h,J) is the log normalisation constant:
  double A = log(pSum);

  // calculate the loglikelihood:
  for(unsigned long long int t = 0; t < nrow; t++) {
    double dataEntry = 0;
    // calculate the loglikelihood of the t'th sample in the data set:
    for (int u = 0; u < d; u++) {
      for(int v = 0; v < d; v++) {
        dataEntry += 0.5 * J(u, v) * data(t, u) * data(t, v);
      }
      dataEntry += h(u) * data(t, u);
    }
    dataEntry -= A;
    // add to the total loglikelihood:
    logLikelihood += dataEntry;
  }

  return(logLikelihood);
}
