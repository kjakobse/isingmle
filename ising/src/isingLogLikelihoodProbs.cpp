#include <Rcpp.h>
using namespace Rcpp;

//' Calculate loglikelihood for data set given distribution
//'
//' \code{isingLogLikelihoodProbs} calculates the value of the loglikelihood function for a given data set in the distribution p.
//'
//' @param data Matrix containing samples from d binary variables in the rows.
//' @param p Numeric vector specifying the distribution of d binary variables.
//' @return \code{isingLogLikelihoodProbs} returns the numeric value of the loglikelihood function.
//' @export
//'
// [[Rcpp::export]]
double isingLogLikelihoodProbs(NumericMatrix data, NumericVector p) {
  int d = data.ncol();
  unsigned long long int nrow = data.nrow();
  double logLikelihood = 0;
  // For each sample look up the corresponding probability in p and add the log to the loglikelihood:
  for(unsigned long long int t = 0; t < nrow; t++) {
    double pEntry = 0;
    for (int u = 0; u < d; u++) {
      if(data(t, u) == 1) {
        pEntry += pow(2, d-u-1);
      }
    }
    logLikelihood += log(p[pEntry]);
  }

  return(logLikelihood);
}
