#include <Rcpp.h>
using namespace Rcpp;


double getXij1(NumericVector &data, IntegerVector &mp,
              IntegerVector &mi, int ri, int ci){
  double v = 0.0;
  for(int k = mp(ci); k < mp(ci+1); k++){
    if(mi(k) == ri) {
      v = data(k);
      break;
    }
  }
  return v;
}


// [[Rcpp::export]]
NumericMatrix sparse_c(SEXP x) {

  S4 mat(x);
  IntegerVector dims = mat.slot("Dim");
  IntegerVector mp = mat.slot("p");
  IntegerVector mi = mat.slot("i");
  NumericVector md = mat.slot("x");
  int nr = dims(0);
  int nc = dims(1);

  int j,i;
  double v;
  NumericMatrix out(nr,nc);
  for(i = 0; i < nr; i++){
    for(j = 0; j < nc; j++) {
      v = getXij1(md, mp, mi, i, j );
      if(v==0) Rprintf(".    ");
      else Rprintf("%4.2f ", v);
      out(i,j) = v;
    }
    Rprintf("\n");
  }

  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
