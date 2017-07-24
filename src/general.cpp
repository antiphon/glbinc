// general stuff
#include "general.h"


// for the sparseMatrix class
double crossprod_std_w_sparse(NumericVector &md, IntegerVector &mi, IntegerVector &mp,
                              NumericVector& center, NumericVector& scale,
                              NumericVector& y, NumericVector &w, int l) {
  double val = 0.0, div = 0.0, xs;
  for(int i = 0; i < y.size(); i++) { // do in two parts
    xs = -center(l);
    val += xs * y(i);
    div += xs * xs * w(i);
  }
  for(int k = mp(l); k < mp(l+1); k++) { // correct for non-zeros
    xs = md(k);
    val += xs * y(mi(k));
    div += (xs * xs - 2 * xs*center(l)) * w(mi(k));
  }
  val /= scale(l);
  div /= scale(l) * scale(l);
  return(val/div);
}

// for the sparseMatrix class. This is slow, try to avoid using.
double getXij(NumericVector &data,
              IntegerVector &mi,
              IntegerVector &mp,
              int ri, int ci){
  double v = 0.0;
  for(int k = mp(ci); k < mp(ci+1); k++){
    if(mi(k) == ri) {
      v = data(k);
      break;
    }
  }
  return v;
}


bool checkConvergence(NumericMatrix& beta, NumericVector& betaold, int l, double eps) {
  int j, converged = 1;
  for(j = 0; j < beta.nrow(); j++) {
    // if(beta(j,l) != 0 & betaold(j) != 0) {
    //   if( std::abs(beta(j,l) - betaold(j)) > eps  ){
    //     converged = 0;
    //     break;
    //   }
    // } else if( beta(j,l) == 0 & betaold(j) != 0){
    //   converged = 0;
    //   break;
    // } else if(beta(j,l)!=0 & betaold(j) == 0) {
    //   converged = 0;
    //   break;
    // }
    if( std::abs(beta(j,l) - betaold(j)) > eps  ){
      converged = 0;
      break;
    }
  }
  if(converged) return true;
  return false;
}

double norm(NumericVector& z) {
  double s=0;
  for(int i=0; i < z.size(); i++) s= s + z(i)*z(i);
  s = sqrt(s);
  return s;
}

double S(double a, double b){
  if(a > b) return a-b;
  if(a < -b) return a+b;
  return 0;
}

// Firm-thresholding operator (from grpreg)
double F(double z, double l1, double l2, double gamma) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (std::abs(z) <= l1) return(0);
  else if (std::abs(z) <= gamma*l1*(1+l2)) return(s*(std::abs(z)-l1)/(1+l2-1/gamma));
  else return(z/(1+l2));
}
// SCAD-modified firm-thresholding operator (from grpreg)
double Fs(double z, double l1, double l2, double gamma) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (std::abs(z) <= l1) return(0);
  else if (std::abs(z) <= (l1*(1+l2)+l1)) return(s*(std::abs(z)-l1)/(1+l2));
  else if (std::abs(z) <= gamma*l1*(1+l2)) return(s*(std::abs(z)-gamma*l1/(gamma-1))/(1-1/(gamma-1)+l2));
  else return(z/(1+l2));
}

