#include <Rcpp.h>
using namespace Rcpp;

//' Updates the distribution of the Ising model and the fitted graph
//'
//' \code{calculateNewPAndEHatnoMTP2} is intended to be used by \code{IsingMLE}. It updates the distribution p for the edges in ePlus when the empirical distributions don't contain zeroes.
//'
//' @param ePlus NumericMatrix containing the edges of G to update p over.
//' @param d Integer containing the number of binary variables.
//' @param e NumericMatrix containing the empirical distributions for the edges in ePlus.
//' @param p NumericVector containing the distribution to be updated.
//' @param eHat NumericMatrix containing the edges in the independence graph of the distribution p.
//' @param capitalJ NumericMatrix containing the canonical parameter J
//' @return \code{calculateNewPAndEHatnoMTP2} returns a list containing the updated distribution its independence graph and the canonical parameter J.
//' @export
//'
// [[Rcpp::export]]
List calculateNewPAndEHatnoMTP2(NumericMatrix ePlus, int d, NumericMatrix e, NumericVector p, NumericMatrix eHat, NumericMatrix capitalJ) {
  int ePlusDim = ePlus.nrow();
  unsigned long long int pLength = p.size();
  for (int t = 0; t < ePlusDim; t++) {
    int i = ePlus(t, 0);
    int j = ePlus(t, 1);
    int iCInternal = d - i;
    int jCInternal = d - j;
    unsigned long long int iMask = 1 << iCInternal;
    unsigned long long int jMask = 1 << jCInternal;

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

    double q11 = e(t, 0) / p11;
    double q10 = e(t, 1) / p10;
    double q01 = e(t, 2) / p01;
    double q00 = e(t, 3) / p00;

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

    unsigned long long int v1 = 0;
    unsigned long long int v2 = iMask;
    unsigned long long int v3 = jMask;
    unsigned long long int v4 = iMask | jMask;
    capitalJ(i-1, j-1) = 0.25 * log(p[v1] * p[v4] / (p[v2] * p[v3]));
    capitalJ(j-1, i-1) = 0.25 * log(p[v1] * p[v4] / (p[v2] * p[v3]));

    if(abs(capitalJ(i-1, j-1)) < 1e-15) {
      eHat(t, 0) = NA_REAL;
      eHat(t, 1) = NA_REAL;
    } else {
      eHat(t, 0) = i;
      eHat(t, 1) = j;
    }
  }

  return List::create(_["p"] = p, _["eHat"] = eHat, _["J"] = capitalJ);
}
