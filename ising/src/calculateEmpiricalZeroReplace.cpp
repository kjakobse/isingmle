#include <Rcpp.h>
using namespace Rcpp;

//' Calculate empirical marginal distributions of pairs in ePlus adding small
//' epsilon to zeroes.
//'
//' \code{CalculateEmpiricalReplaceZeroes} calculates the empirical distribution
//' of pairs of binary variables specified by ePlus from the sample first and
//' second moments and adds a small positive constant if a probability is zero.
//'
//' ePlus should contain the variable pairs for which the empirical distribution
//' is to be calculated. They should be coded as integers from 1 up to the
//' number of binary variables.
//' The sample moments should be ordered according to this coding, so the sample
//' mean of variable 1 should have index 1 in xBar (using 1-indexing).
//' The value added to the zeroes is subsequently subtracted from the positive
//' values so they still sum to 1.
//'
//' @param ePlus Matrix containing pairs of variables in the rows.
//' @param m Matrix containing the sample second moment.
//' @param xBar Numeric vector containing the sample first moment.
//' @param ReplaceValue Double to replace zeroes with. Must be >0.
//' @return \code{CalculateEmpirical} returns a matrix with the calculated
//' empirical distributions.
//' @author Kim Daniel Jakobsen
//' @examples
//' 1+1
//'
//' @export
//'
// [[Rcpp::export]]

NumericMatrix calculateEmpiricalReplaceZeroes(NumericMatrix ePlus,
                                              NumericMatrix m,
                                              NumericVector xBar,
                                              double ReplaceValue) {
  // determine the number of pairs to calculate the empirical distribution of:
  int ePlusDim = ePlus.nrow();
  // initiate the integers counting how many observations in the pair-margin
  // are zero and positive:
  int zerocounter = 0;
  int poscounter = 0;
  NumericMatrix empirical(ePlusDim, 4);
  for (int t = 0; t < ePlusDim; t++) {
    // for each row in ePlus, extract the variable pair, and calculate the
    // empirical distribution from the corresponding entries in the sample
    // moments:
    int i = ePlus(t, 0) - 1;
    int j = ePlus(t, 1) - 1;
    empirical(t, 0) = (1 + xBar[i] + xBar[j] + m(i, j)) / 4;
    empirical(t, 1) = (1 + xBar[i] - xBar[j] - m(i, j)) / 4;
    empirical(t, 2) = (1 - xBar[i] + xBar[j] - m(i, j)) / 4;
    empirical(t, 3) = (1 - xBar[i] - xBar[j] + m(i, j)) / 4;
    // Check that none of the probabilities are negative:
    if (empirical(t, 0) < -1e-15 ||
        empirical(t, 1) < -1e-15 ||
        empirical(t, 2) < -1e-15 ||
        empirical(t, 3) < -1e-15) {
      NumericMatrix error = NumericMatrix(0, 0);
      return(error);
    }
    // check which probabilities are zero, if so replace with ReplaceValue and
    // increment the zero-counter, otherwise increment positive-counter:
    if (empirical(t, 0) <= 0) {
      empirical(t, 0) = ReplaceValue;
      zerocounter++;
    } else if(empirical(t, 0) > ReplaceValue) {
      poscounter++;
    }
    if (empirical(t, 1) <= 0) {
      empirical(t, 1) = ReplaceValue;
      zerocounter++;
    } else if(empirical(t, 1) > ReplaceValue) {
      poscounter++;
    }
    if (empirical(t, 2) <= 0) {
      empirical(t, 2) = ReplaceValue;
      zerocounter++;
    } else if(empirical(t, 2) > ReplaceValue) {
      poscounter++;
    }
    if (empirical(t, 3) <= 0) {
      empirical(t, 3) = ReplaceValue;
      zerocounter++;
    } else if(empirical(t, 3) > ReplaceValue) {
      poscounter++;
    }
    // check which probabilities are positive, and subtract the number of
    // ReplaceValue's which are added to the zeroes:
    if (empirical(t, 0) > ReplaceValue) {
      empirical(t, 0) -= ReplaceValue * zerocounter / poscounter;
    }
    if (empirical(t, 1) > ReplaceValue) {
      empirical(t, 1) -= ReplaceValue * zerocounter / poscounter;
    }
    if (empirical(t, 2) > ReplaceValue) {
      empirical(t, 2) -= ReplaceValue * zerocounter / poscounter;
    }
    if (empirical(t, 3) > ReplaceValue) {
      empirical(t, 3) -= ReplaceValue * zerocounter / poscounter;
    }
    // make sure the empirical probabilities are normalised:
    double empiricalsum = empirical(t, 0) +
                          empirical(t, 1) +
                          empirical(t, 2) +
                          empirical(t, 3);
    if (empiricalsum != 1) {
      empirical(t, 0) /= empiricalsum;
      empirical(t, 1) /= empiricalsum;
      empirical(t, 2) /= empiricalsum;
      empirical(t, 3) /= empiricalsum;
    }
  }
  return(empirical);
}
