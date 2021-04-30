#include <Rcpp.h>
using namespace Rcpp;

//' calculateNewPAndEHatnoMTP2Boundary <to be documented>
//'
//' @param ePlus <to be documented>
//' @param d <to be documented>
//' @param e <to be documented>
//' @param p <to be documented>
//' @param eHat <to be documented>
//' @return List <to be documented>
//' @export
//'
// [[Rcpp::export]]
List calculateNewPAndEHatnoMTP2Boundary(NumericMatrix ePlus, int d, NumericMatrix e, NumericVector p, NumericMatrix eHat) {
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

    double q11;
    double q10;
    double q01;
    double q00;

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
    eHat(t, 0) = i;
    eHat(t, 1) = j;

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
