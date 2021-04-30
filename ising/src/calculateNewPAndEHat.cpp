#include <Rcpp.h>
using namespace Rcpp;

//' calculateNewPAndEHat <to be documented>
//'
//' @param ePlus <to be documented>
//' @param d <to be documented>
//' @param e <to be documented>
//' @param p <to be documented>
//' @param eHat <to be documented>
//' @param capitalJ <to be documented>
//' @param iter <to be documented>
//' @return List <to be documented>
//' @export
//'
// [[Rcpp::export]]
List calculateNewPAndEHat(NumericMatrix ePlus, int d, NumericMatrix e, NumericVector p, NumericMatrix eHat, NumericMatrix capitalJ, int iter) {
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

    double delta = 0.25 * log(q11 * q00 / (q10 * q01));

    unsigned long long int v1 = 0;
    unsigned long long int v2 = iMask;
    unsigned long long int v3 = jMask;
    unsigned long long int v4 = iMask | jMask;
    if(iter == 0) {
      capitalJ(i-1, j-1) = 0.25 * log(p[v1] * p[v4] / (p[v2] * p[v3]));
      capitalJ(j-1, i-1) = 0.25 * log(p[v1] * p[v4] / (p[v2] * p[v3]));
    }

    double v00;
    double v10;
    double v01;
    double v11;
    if (delta + capitalJ(i-1, j-1) > 0) { //changed >= to >
      v11 = q11;
      v10 = q10;
      v01 = q01;
      v00 = q00;
      eHat(t, 0) = i;
      eHat(t, 1) = j;
    } else {
      double capitalC = exp(-4 * capitalJ(i-1, j-1)) * (p11 * p00) / (p10 * p01);
      double a = 1 - capitalC;
      double b = e(t, 0) + e(t, 3) + capitalC * (e(t, 1) + e(t, 2));
      double c = e(t, 0) * e(t, 3) - capitalC * (e(t, 1) * e(t, 2));
      double lambda;
      if(a == 0) {
        lambda = -c / b;
      } else {
        lambda = (-b + sqrt(pow(b,2) - 4 * a * c)) / (2 * a);
      }
      v11 = (e(t, 0) + lambda) / p11;
      v10 = (e(t, 1) - lambda) / p10;
      v01 = (e(t, 2) - lambda) / p01;
      v00 = (e(t, 3) + lambda) / p00;
      eHat(t, 0) = NA_REAL;
      eHat(t, 1) = NA_REAL;
    }

    for (unsigned long long int u = 0; u < pLength; u++) {
      if (u & iMask) {
        if (u & jMask) {
          p[u] = p[u] * v11;
        } else {
          p[u] = p[u] * v10;
        }
      } else {
        if (u & jMask) {
          p[u] = p[u] * v01;
        } else {
          p[u] = p[u] * v00;
        }
      }
    }

    capitalJ(i-1, j-1) = 0.25 * log(p[v1] * p[v4] / (p[v2] * p[v3]));
    capitalJ(j-1, i-1) = 0.25 * log(p[v1] * p[v4] / (p[v2] * p[v3]));
  }

  return List::create(_["p"] = p, _["eHat"] = eHat, _["J"] = capitalJ);
}
