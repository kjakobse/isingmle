#include <Rcpp.h>
using namespace Rcpp;

//' calculateNewPAndEHatBoundary <to be documented>
//'
//' @param ePlus <to be documented>
//' @param d <to be documented>
//' @param e <to be documented>
//' @param p <to be documented>
//' @param eHat <to be documented>
//' @return List <to be documented>
//' @export
//
// [[Rcpp::export]]
List calculateNewPAndEHatBoundary(NumericMatrix ePlus, int d, NumericMatrix e, NumericVector p, NumericMatrix eHat) {
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

    if (abs(q11) < 1e-15 || abs(q10) < 1e-15 || abs(q01) < 1e-15 || abs(q00) < 1e-15) {
      eHat(t, 0) = i;
      eHat(t, 1) = j;
    } else {
      double delta = 0.25 * log(q11 * q00 / (q10 * q01));

      unsigned long long int v1 = 0;
      unsigned long long int v2 = iMask;
      unsigned long long int v3 = jMask;
      unsigned long long int v4 = iMask | jMask;

      double capitalJ = 0;
      unsigned long long int wLength = pow(2, d-2);
      unsigned long long int index = 0;
      for(unsigned long long int w = 0; w < wLength; w++) {
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
        index++;
      }

      if (delta + capitalJ > 0) { //changed >= to >
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
          lambda = (-b + sqrt(pow(b,2) - 4 * a * c)) / (2 * a);
        }
        q11 = (e(t, 0) + lambda) / p11;
        q10 = (e(t, 1) - lambda) / p10;
        q01 = (e(t, 2) - lambda) / p01;
        q00 = (e(t, 3) + lambda) / p00;
        eHat(t, 0) = NA_REAL;
        eHat(t, 1) = NA_REAL;
      }
    }

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
