#include <Rcpp.h>
using namespace Rcpp;

//' Updates the distribution of the Ising model on the boundary
//'
//' \code{calculateNewPAndEHatnoMTP2Boundary} is intended to be used by \code{IsingMLE}. It updates the distribution p for the edges in ePlus when some of the empirical distributions contain zeroes.
//'
//' @param ePlus Matrix containing the edges of G to update p over.
//' @param d Integer containing the number of binary variables.
//' @param e Matrix containing the empirical distributions for the edges in ePlus.
//' @param p Numeric vector containing the distribution to be updated.
//' @return \code{calculateNewPAndEHatnoMTP2Boundary} returns a list containing the updated distribution.
//' @export
//'
// [[Rcpp::export]]
List calculateNewPAndEHatnoMTP2Boundary(NumericMatrix ePlus,
                                        int d,
                                        NumericMatrix e,
                                        NumericVector p) {
  int ePlusDim = ePlus.nrow();
  unsigned long long int pLength = p.size();
  // for loop over the edges in ePlus:
  for (int t = 0; t < ePlusDim; t++) {
    int i = ePlus(t, 0);
    int j = ePlus(t, 1);
    // masks for the variables i and j. << is the bitwise and operator and shifts the bits to the left:
    int iCInternal = d - i;
    int jCInternal = d - j;
    unsigned long long int iMask = 1 << iCInternal;
    unsigned long long int jMask = 1 << jCInternal;

    // calculate the marginal distribution of variables i and j:
    double p00 = 0.0;
    double p01 = 0.0;
    double p10 = 0.0;
    double p11 = 0.0;
    for (unsigned long long int u = 0; u < pLength; u++) {
      if (u & iMask) {
        if (u & jMask) {
          p11 += p[u];
        } else {
          p10 += p[u];
        }
      } else {
        if (u & jMask) {
          p01 += p[u];
        } else {
          p00 += p[u];
        }
      }
    }

    double q11;
    double q10;
    double q01;
    double q00;

    // calculate the factors with which to update the probabilities in p. If the empirical distribution is 0 the corresponding q is set equal to 0:
    if(abs(e(t, 0) < 1e-15) || abs(p11) < 1e-15) {
      q11 = 0;
    } else {
      q11 = e(t, 0) / p11;
    }
    if(abs(e(t, 1) < 1e-15) || abs(p10) < 1e-15) {
      q10 = 0;
    } else {
      q10 = e(t, 1) / p10;
    }
    if(abs(e(t, 2) < 1e-15) || abs(p01) < 1e-15) {
      q01 = 0;
    } else {
      q01 = e(t, 2) / p01;
    }
    if(abs(e(t, 3) < 1e-15) || abs(p00) < 1e-15) {
      q00 = 0;
    } else {
      q00 = e(t, 3) / p00;
    }

    // each entry in p is updated by checking the value of the i'th and j'th variable using the masks and multiplying with the appropriate q:
    for (unsigned long long int u = 0; u < pLength; u++) {
      if (u & iMask) {
        if (u & jMask) {
          p[u] = p[u] * q11;
        } else {
          p[u] = p[u] * q10;
        }
      } else {
        if (u & jMask) {
          p[u] = p[u] * q01;
        } else {
          p[u] = p[u] * q00;
        }
      }
    }
  }
  return List::create(_["p"] = p);
}
