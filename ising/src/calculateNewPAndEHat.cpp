#include <Rcpp.h>
using namespace Rcpp;

//' Updates the distribution of the Ising model with MTP2 constraints and the fitted graph
//'
//' \code{calculateNewPAndEHat} is intended to be used by \code{IsingMLEmtp2}. It updates the distribution p for the edges in ePlus when the empirical distributions don't contain zeroes.
//'
//' \code{calculateNewPAndEHat} doesn't use the value of eHat. It fills the matrix eHat with the edges of the independence graph after the update.\cr
//' iter is used to indicate if the provided J should be used or not. If iter = 0 J is calculated and otherwise the provided J is used.
//'
//' @param ePlus NumericMatrix containing the edges of G to update p over.
//' @param d Integer containing the number of binary variables.
//' @param e NumericMatrix containing the empirical distributions for the edges in ePlus.
//' @param p NumericVector containing the distribution to be updated.
//' @param eHat NumericMatrix containing the edges in the independence graph of the distribution p.
//' @param capitalJ NumericMatrix containing the canonical parameter J.
//' @param iter Integer specifying the iteration of the algorithm.
//' @return \code{calculateNewPAndEHat} returns a list containing the updated distribution, its independence graph, and the canonical parameter J.
//' @export
//'
// [[Rcpp::export]]
List calculateNewPAndEHat(NumericMatrix ePlus, int d, NumericMatrix e, NumericVector p, NumericMatrix eHat, NumericMatrix capitalJ, int iter) {
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

// calculate the factors with which to update the probabilities in p if delta+capitalJ>0:
    double q11 = e(t, 0) / p11;
    double q10 = e(t, 1) / p10;
    double q01 = e(t, 2) / p01;
    double q00 = e(t, 3) / p00;

// calculate delta:
    double delta = 0.25 * log(q11 * q00 / (q10 * q01));

// Define indices used to calculate capitalJ:
    unsigned long long int v1 = 0;
    unsigned long long int v2 = iMask;
    unsigned long long int v3 = jMask;
    unsigned long long int v4 = iMask | jMask;
// If iter = 0, calculate J:
    if(iter == 0) {
      capitalJ(i-1, j-1) = 0.25 * log(p[v1] * p[v4] / (p[v2] * p[v3]));
      capitalJ(j-1, i-1) = 0.25 * log(p[v1] * p[v4] / (p[v2] * p[v3]));
    }

// If delta+J_ij>0, leave the q's unchanged, otherwise shift them so J_ij = 0:
    double v00;
    double v10;
    double v01;
    double v11;
    if (delta + capitalJ(i-1, j-1) > 0) {
      v11 = q11;
      v10 = q10;
      v01 = q01;
      v00 = q00;
      // when the updated J_ij is positive the edge ij is in the independence graph:
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
        lambda = (-b + sqrt(pow(b,2) - 4 * a * c)) / (2 * a); // the shift is the root of a quadratic function.
      }
      v11 = (e(t, 0) + lambda) / p11;
      v10 = (e(t, 1) - lambda) / p10;
      v01 = (e(t, 2) - lambda) / p01;
      v00 = (e(t, 3) + lambda) / p00;
      // when the updated J_ij is zero the edge ij is not in the independence graph:
      eHat(t, 0) = NA_REAL;
      eHat(t, 1) = NA_REAL;
    }

// each entry in p is updated by checking the value of the i'th and j'th variable using the masks and multiplying with the appropriate q:
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

// The value of the canonical parameter J is calculated for the updated distribution:
    capitalJ(i-1, j-1) = 0.25 * log(p[v1] * p[v4] / (p[v2] * p[v3]));
    capitalJ(j-1, i-1) = 0.25 * log(p[v1] * p[v4] / (p[v2] * p[v3]));
  }

  return List::create(_["p"] = p, _["eHat"] = eHat, _["J"] = capitalJ);
}
