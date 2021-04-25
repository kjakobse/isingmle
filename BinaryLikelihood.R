library(Rcpp)

cppFunction('
  double binaryLogLikelihood(NumericVector intData, int d) {
    unsigned long long int length = intData.size();
    double logLikelihood = 0;
    NumericVector counts(pow(2,d));
    
    for(unsigned long long int t = 0; t < length; t++) {
      counts[intData[t]]++;
    }
    
    for(unsigned long long int i = 0; i < pow(2, d); i++) {
      if(counts[i] > 0) {
        logLikelihood += counts[i] * log(counts[i] / length);
      }
    }
    
    return(logLikelihood);
  }
')
