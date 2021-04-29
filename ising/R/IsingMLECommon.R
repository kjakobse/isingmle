library(Rcpp)

# cppFunction('
#   NumericVector calculatexBar(NumericMatrix data, unsigned long long int length, int d) {
#     NumericVector xBar(d);
#
#     for (int t = 0; t < d; t++) {
#       double dEntry = 0;
#
#       for (unsigned long long int u = 0; u < length; u++) {
#         dEntry += data(u, t);
#       }
#       xBar[t] = dEntry / length;
#     }
#
#     return (xBar);
#   }
# ')
#
# cppFunction('
#   NumericMatrix calculateM(NumericMatrix data, unsigned long long int length, int d) {
#     NumericMatrix matrixM(d, d);
#
#     for (int t = 0; t < d; t++) {
#       for (int u = 0; u < d; u++) {
#         double entryValue = 0;
#         for (unsigned long long int v = 0; v < length; v++) {
#           //entryValue += data(v, t) * data(v, u);
#           if (data(v, t) == data(v, u)) {
#             entryValue++;
#           } else {
#             entryValue--;
#           }
#         }
#         matrixM(t, u) = entryValue / length;
#       }
#     }
#
#     return(matrixM);
#   }
# ')
#
# cppFunction('
#   NumericVector createP(NumericVector mu, int d) {
#     unsigned long long int pLength = pow(2, d);
#     NumericVector p(pLength);
#     double normalisationConstant = pow(2, -d);
#     for (unsigned long long int t = 0; t < pLength; t++) {
#       double pEntry = 1;
#       for (int u = 0; u < d; u++) {
#         unsigned long long int tMask = 1 << u;
#         int observationValue = -1;
#         if (t & tMask) {
#           observationValue = 1;
#         }
#         double value = 1 + observationValue * mu[d - u - 1];
#         pEntry *= value;
#       }
#       p[t] = normalisationConstant * pEntry;
#     }
#
#     return(p);
#   }
# ')
#
# cppFunction('
#   NumericVector calculateNewMu(NumericVector p, int d) {
#     NumericVector mu(d);
#     unsigned long long int pLength = p.size();
#
#     for (int t = 0; t < d; t++) {
#       unsigned long long int uMask = 1 << (d - t - 1);
#
#       double dEntry = 0;
#       for (unsigned long long int u = 0; u < pLength; u++) {
#         if (u & uMask) {
#           dEntry += p[u]; // Observation is +1 add p[u]
#         } else {
#           dEntry -= p[u]; // observation is -1 subtract p[u]
#         }
#       }
#       mu[t] = dEntry;
#     }
#
#     return(mu);
#   }
# ')
#
# cppFunction('
#   NumericMatrix calculateNewXi(NumericVector p, int d) {
#     NumericMatrix xi(d, d);
#     unsigned long long int pLength = p.size();
#
#     for (int t = 0; t < d; t++) {
#       unsigned long long int tMask = 1 << (d - t - 1);
#       for (int u = t; u < d; u++) {
#         unsigned long long int uMask = 1 << (d - u - 1);
#
#         double entryValue = 0;
#         for (unsigned long long int v = 0; v < pLength; v++) {
#           if (((bool)(v & tMask)) == ((bool)(v & uMask))) {
#             entryValue += p[v]; // observations have same value
#           } else {
#             entryValue -= p[v]; // observations have opposite value
#           }
#         }
#         xi(t, u) = entryValue;
#         if (t != u) { // save symmetric value
#           xi(u, t) = entryValue;
#         }
#       }
#     }
#
#     return(xi);
#   }
# ')
#
# cppFunction('
#   NumericMatrix createEPlus(NumericMatrix e, NumericMatrix m, NumericVector xBar) {
#     int nRows = e.nrow();
#     int nCols = e.ncol();
#
#     NumericVector includeRows(nRows);
#     int nIncludedRows = 0;
#     for (int i = 0; i < nRows; i++) {
#       if (m(e(i, 0) - 1, e(i, 1) - 1) > (xBar[e(i, 0) - 1] * xBar[e(i, 1) - 1])) {
#         includeRows[i] = 1;
#         nIncludedRows++;
#       }
#     }
#     NumericMatrix ePlus(nIncludedRows, nCols);
#     int k = 0;
#     for (int i = 0; i < nRows; i++) {
#       if (includeRows[i] == 1) {
#         for (int j = 0; j < nCols; j++) {
#           ePlus(k, j) = e(i, j);
#         }
#         k++;
#       }
#     }
#
#     return(ePlus);
#   }
# ')
#
# cppFunction('
#   NumericVector calculateCondition(NumericMatrix eHatOmitNA, NumericMatrix m, NumericMatrix xi) {
#     int nRows = eHatOmitNA.nrow();
#
#     NumericVector condition(nRows);
#     for (int i = 0; i < nRows; i++) {
#       condition[i] = abs(m(eHatOmitNA(i, 0) - 1, eHatOmitNA(i, 1) - 1) - xi(eHatOmitNA(i, 0) - 1, eHatOmitNA(i, 1) - 1));
#     }
#
#     return(condition);
#   }
# ')
#
# cppFunction('
#   NumericVector calculateCondition2(NumericMatrix e, NumericMatrix m, NumericMatrix xi) {
#     int nRows = e.nrow();
#
#     NumericVector condition2(nRows);
#     for (int i = 0; i < nRows; i++) {
#       condition2[i] = m(e(i, 0) - 1, e(i, 1) - 1) - xi(e(i, 0) - 1, e(i, 1) - 1);
#     }
#
#     return(condition2);
#   }
# ')
#
# cppFunction('
#   List CalculateEmpirical(NumericMatrix ePlus, NumericMatrix m, NumericVector xBar) {
#     int ePlusDim = ePlus.nrow();
#     int containsZeroes = 0;
#
#     NumericMatrix empirical(ePlusDim, 4);
#     for (int t = 0; t < ePlusDim; t++) {
#       int i = ePlus(t, 0) - 1;
#       int j = ePlus(t, 1) - 1;
#       empirical(t, 0) = (1 + xBar[i] + xBar[j] + m(i, j)) / 4;
#       empirical(t, 1) = (1 + xBar[i] - xBar[j] - m(i, j)) / 4;
#       empirical(t, 2) = (1 - xBar[i] + xBar[j] - m(i, j)) / 4;
#       empirical(t, 3) = (1 - xBar[i] - xBar[j] + m(i, j)) / 4;
#       if (abs(empirical(t, 0)) <= 1e-15 || abs(empirical(t, 1)) <= 1e-15 || abs(empirical(t, 2)) <= 1e-15 || abs(empirical(t, 3)) <= 1e-15) {
#         containsZeroes = 1;
#       }
#       if (empirical(t, 0) < -1e-15 || empirical(t, 1) < -1e-15 || empirical(t, 2) < -1e-15 || empirical(t, 3) < -1e-15) {
#         NumericMatrix error = NumericMatrix(0, 0);
#         return List::create(_["empirical"] = error);
#       }
#     }
#     return List::create(_["empirical"] = empirical, _["containsZeroes"] = containsZeroes);
#   }
# ')
#
# cppFunction('
#   NumericVector calculateH(int d, NumericVector p) {
#     NumericVector h(d);
#
#     unsigned long long int p2powd = 1 << d;
#     for (int i = 0; i < d; i++) {
#       unsigned long long int p2powDMinusI = 1 << (d-i-1);
#       h[i] = log(p[p2powDMinusI] * p[p2powd - 1] / (p[p2powd - p2powDMinusI - 1] * p[0])) / 4;
#     }
#
#     return(h);
#   }
# ')
