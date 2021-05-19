#include <Rcpp.h>
using namespace Rcpp;

//' Updates the distribution of the Ising model with MTP2 constraints on the boundary
//'
//' \code{calculateNewPAndEHatBoundary} is intended to be used by \code{IsingMLEmtp2}. It updates the distribution p for the edges in ePlus when some of the empirical distributions contain zeroes.
//'
//' @param ePlus Matrix containing the edges of G to update p over.
//' @param d Integer containing the number of binary variables.
//' @param e Matrix containing the empirical distributions for the edges in ePlus.
//' @param p Numeric vector containing the distribution to be updated.
//' @param eHat Matrix containing the edges in the independence graph of the distribution p.
//' @return \code{calculateNewPAndEHatBoundary} returns a list containing the updated distribution and its independence graph.
//' @export
//'
// [[Rcpp::export]]
List calculateNewPAndEHatBoundary(NumericMatrix ePlus, int d, NumericMatrix e, NumericVector p, NumericMatrix eHat) {
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
    if(abs(e(t, 2)) < 1e-15 || abs(p01) < 1e-15) {
      q01 = 0;
    } else {
      q01 = e(t, 2) / p01;
    }
    if(abs(e(t, 3)) < 1e-15 || abs(p00) < 1e-15) {
      q00 = 0;
    } else {
      q00 = e(t, 3) / p00;
    }
// If one of the q's is zero they remain unchanged, otherwise we attempt to calculate J_ij:
    if (abs(q11) < 1e-15 || abs(q10) < 1e-15 || abs(q01) < 1e-15 || abs(q00) < 1e-15) {
      eHat(t, 0) = i;
      eHat(t, 1) = j;
    } else {
      // calculate delta:
      double delta = 0.25 * log(q11 * q00 / (q10 * q01));
      // Define indices used to calculate J_ij:
      unsigned long long int v1 = 0;
      unsigned long long int v2 = iMask;
      unsigned long long int v3 = jMask;
      unsigned long long int v4 = iMask | jMask;
      // calculate J_ij by running through all combinations of the variables other than i and j until one has all probabilities positive:
      double capitalJ = 0;
      unsigned long long int indexLength = pow(2, d-2);
      for(unsigned long long int index = 0; index < indexLength; index++) {
        if(jMask & index) {
          index += jMask;
        }
        if(iMask & index) {
          index += iMask;
        }
        if(jMask & index) {
          index += jMask;
        }
        if(p[v1 + index] > 1e-15 && p[v2 + index] > 1e-15 && p[v3 + index] > 1e-15 && p[v4 + index] > 1e-15) {
          capitalJ = 0.25 * log(p[v1 + index] * p[v4 + index] / (p[v2 + index] * p[v3 + index]));
          break;
        }
        // If no indices allow us to calculate J_ij return an error:
        if(index == (indexLength - 1)) {
          return List::create(_["p"] = NA_REAL);
        }
      }
      // If delta+J_ij>0, leave the q's unchanged, otherwise shift them so J_ij = 0:
      if (delta + capitalJ > 0) {
        // when the updated J_ij is positive the edge ij is in the independence graph:
        eHat(t, 0) = i;
        eHat(t, 1) = j;
      } else {
        double capitalC = exp(-4 * capitalJ) * (p11 * p00) / (p10 * p01);
        double a = 1 - capitalC;
        double b = e(t, 0) + e(t, 3) + capitalC * (e(t, 1) + e(t, 2));
        double c = e(t, 0) * e(t, 3) - capitalC * (e(t, 1) * e(t, 2));
        double lambda;
        if(a == 0) {
          lambda = -c / b;
        } else {
          lambda = (-b + sqrt(pow(b,2) - 4 * a * c)) / (2 * a); // the shift is the root of a quadratic function.
        }
        q11 = (e(t, 0) + lambda) / p11;
        q10 = (e(t, 1) - lambda) / p10;
        q01 = (e(t, 2) - lambda) / p01;
        q00 = (e(t, 3) + lambda) / p00;
        // when the updated J_ij is zero the edge ij is not in the independence graph:
        eHat(t, 0) = NA_REAL;
        eHat(t, 1) = NA_REAL;
      }
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

  return List::create(_["p"] = p, _["eHat"] = eHat);
}
