#include <Rcpp.h>
using namespace Rcpp;

double crossprod_std_w_sparse(NumericVector &md, IntegerVector &mi, IntegerVector &mp,
                              NumericVector& center, NumericVector& scale,
                              NumericVector& y, NumericVector &w, int l);
double getXij(NumericVector &data, IntegerVector &mi, IntegerVector &mp, int ri, int ci);
bool checkConvergence(NumericMatrix& beta, NumericVector& betaold, int l, double eps) ;
double norm(NumericVector& z) ;
double S(double a, double b);
