#include <Rcpp.h>
using namespace Rcpp;

//' Calculate the log-likelihood for the MLE among all binary distributions.
//'
//' \code{binaryLogLikelihood} calculates the log-likelihood for a data set
//' from d binary variables in the maximum likelihood estimator among all
//' binary distributions.
//'
//' The data must be given as a vector of integers, where each entry
//' corresponds to an observation from the d binary variables with the value of
//'  its binary representation equal to the value of the integer.
//'
//' @param intData NumericVector containing the data encoded as integers.
//' @param d Integer specifying the number of binary variables.
//' @return \code{binaryLogLikelihood} returns a numeic value with the
//' log-likelihood.
//' @export
//'
// [[Rcpp::export]]
double binaryLogLikelihood(NumericVector intData,
                           int d) {
  unsigned long long int length = intData.size();
  double logLikelihood = 0;
  NumericVector counts(pow(2, d));

// count the number of each outcome observed:
  for(unsigned long long int t = 0; t < length; t++) {
    counts[intData[t]]++;
  }

// Calculate the log-likelihood for the model with probabilities equal to
//the sample proportions:
  for(unsigned long long int i = 0; i < pow(2, d); i++) {
    if(counts[i] > 0) {
      logLikelihood += counts[i] * log(counts[i] / length);
    }
  }

  return(logLikelihood);
}
