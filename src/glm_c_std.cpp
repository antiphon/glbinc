#include <Rcpp.h>
#include <Rmath.h>
#include "general.h"
using namespace Rcpp;


// update per one coefficient v'x_j/x_j'Wx_j
double crossprod_std_w(NumericMatrix &x, NumericVector& center, NumericVector& scale,
                     NumericVector& y, NumericVector &w, int l) {
  double val = 0, div= 0, xs;
  for (int i=0; i< y.size(); i++){
    xs = (x(i,l) - center(l))/scale(l);
    val += xs * y(i);
    div += xs * xs * w(i);
  }
  return(val/div);
}



// [[Rcpp::export]]
List glm_binom_std_c(NumericMatrix X,
                     NumericVector y,
                     NumericVector offset,
                     NumericVector center,
                     NumericVector scale,
                     bool add_intercept,
                     int verb,
                     double eps,
                     int maxiter,
                     double w_limit = 0.25,
                     double stability_threshold = 1e-3
) {

  int n = y.size();
  int p = X.ncol();

  int i, j, it;
  double pi, pip, shift, si;
  double unpenalised_shifts;
  double sum_w, sum_v;

  //
  double a0, b0;
  NumericVector beta(p);
  NumericVector betaold(p);
  NumericVector eta(n);
  NumericVector r(n);
  NumericVector v(n);
  NumericVector w(n);

  // initial
  double ybar = mean(y);
  a0 = log(ybar/(1.0-ybar));
  for(i=0; i < n; i++){
    eta(i) = a0 + offset(i);
  }
  it = 0;

  double threshold = stability_threshold; // for numeric stability

  // main loop
  while(it < maxiter) {
    it++;
    if(verb) Rprintf("        \r[it:%i diff %f]", it, unpenalised_shifts);
    // quadratic approximation of the log likelihood
    sum_w = 0.0;
    sum_v = 0.0;
    unpenalised_shifts = 0.0;
    // update quadratic approx
    for(i=0; i < n; i++) {
      if(eta(i) > 10) {
        pi = 1;
        pip = 1 - threshold;
      }else if(eta(i) < -10){
        pi = 0;
        pip = threshold;
      } else{
        pi = 1.0/(1.0 + exp(-eta(i)));
        pip = pi;
      }
      w(i) = fmax(w_limit, pip * (1.0-pip)); // hessian // Unstable as free
      v(i) = y(i) - pi; // gradient
      r(i) = v(i) / w(i); // residual
      sum_w += w(i);
      sum_v += v(i);
    }
    // update intercept if needed
    if(add_intercept){
      shift = sum_v/sum_w;
      b0 = shift + a0;
      for(i=0; i < n; i++) {
        eta(i) += shift;
      }
      unpenalised_shifts += std::abs(shift);
    }
    // update covariate coefficients
    for(j=0; j < p; j++) {
      shift = crossprod_std_w(X, center, scale, v, w, j);
      unpenalised_shifts += std::abs(shift);
      beta(j) = betaold(j) + shift;
      for(i=0; i < n; i++) {
        si = shift * (X(i,j)-center(j))/scale(j);
        eta(i) += si;
      }
    }

    // check convergence
    if(unpenalised_shifts < eps){
      break;
    }

    //
    a0 = b0;
    for(j = 0; j < p; j++) betaold(j) = beta(j);
  }// end of main loop
  if(verb) Rprintf("\n");
  // done. Update w and r
  for(i=0; i < n; i++) {
    if(eta(i) > 10) {
      pi = 1;
      pip = 1 - threshold;
    }else if(eta(i) < -10){
      pi = 0;
      pip = threshold;
    } else{
      pi = 1.0/(1.0 + exp(-eta(i)));
      pip = pi;
    }
    w(i) = pip * (1.0-pip); // hessian // Unstable as free
    v(i) = y(i) - pi; // gradient
    r(i) = v(i) / w(i); // residual
  }
  //
  return List::create(Named("beta0")=b0, Named("beta")=beta, Named("r")=r, Named("w")=w);
}




// [[Rcpp::export]]
List glm_binom_std_sparse_c(SEXP X,
                            NumericVector y,
                            NumericVector offset,
                            NumericVector center,
                            NumericVector scale,
                            bool add_intercept,
                            int verb,
                            double eps,
                            int maxiter,
                            double w_limit = 0.25,
                            double stability_threshold = 1e-3
) {
  S4 mat(X);
  IntegerVector dims = mat.slot("Dim");
  IntegerVector mp = mat.slot("p");
  IntegerVector mi = mat.slot("i");
  NumericVector md = mat.slot("x");
  int nc = dims(1);
  int n = y.size();
  int p = nc;

  int i, j, k, it, it2;
  double pi, pip, shift, si;
  double unpenalised_shifts;
  double sum_w, sum_v, xs;

  //
  double a0, b0;
  NumericVector beta(p);
  NumericVector betaold(p);
  NumericVector eta(n);
  NumericVector r(n);
  NumericVector v(n);
  NumericVector w(n);

  // initial
  double ybar = mean(y);
  a0 = log(ybar/(1.0-ybar));
  for(i=0; i < n; i++){
    eta(i) = a0 + offset(i);
  }
  it = 0;
  it2 = 0;
  double threshold = stability_threshold; // for numeric stability

  // main loop
  while(it < maxiter) {
    it++;
    if(verb) Rprintf("        \r[it:%i diff %f]", it, unpenalised_shifts);
    // quadratic approximation of the log likelihood
    sum_w = 0.0;
    sum_v = 0.0;
    unpenalised_shifts = 0.0;
    // update quadratic approx
    for(i=0; i < n; i++) {
      if(eta(i) > 10) {
        pi = 1;
        pip = 1 - threshold;
      }else if(eta(i) < -10){
        pi = 0;
        pip = threshold;
      } else{
        pi = 1.0/(1.0 + exp(-eta(i)));
        pip = pi;
      }
      w(i) = fmax(w_limit, pip * (1.0-pip) ) ; // hessian // Unstable as free
      v(i) = y(i) - pi; // gradient
      r(i) = v(i) / w(i); // residual
      sum_w += w(i);
      sum_v += v(i);
    }
    // update intercept if needed
    if(add_intercept){
      shift = sum_v/sum_w;
      b0 = shift + a0;
      for(i=0; i < n; i++) {
        eta(i) += shift;
      }
      unpenalised_shifts += std::abs(shift);
    }
    // update non-intercept coefficients
    for(j=0; j < p; j++) {
      shift = crossprod_std_w_sparse(md, mi, mp, center, scale, v, w, j);
      unpenalised_shifts += std::abs(shift);
      beta(j) = betaold(j) + shift;
      for(i=0; i < n; i++) {
        xs =-center(j)/scale(j);
        si = shift * xs;
        eta(i) += si;
      }
      for(k = mp(j); k < mp(j+1); k++) {
        xs = md(k)/scale(j);
        si = shift * xs;
        eta(mi(k)) += si;
      }

    }
    // check convergence
    if(unpenalised_shifts < eps){
      break;
    }
    //
    a0 = b0;
    for(j = 0; j < p; j++) betaold(j) = beta(j);
  }// end of main loop
  if(verb) Rprintf(" ok\n");
  // done. Update w and r
  for(i=0; i < n; i++) {
    if(eta(i) > 10) {
      pi = 1;
      pip = 1 - threshold;
    }else if(eta(i) < -10){
      pi = 0;
      pip = threshold;
    } else{
      pi = 1.0/(1.0 + exp(-eta(i)));
      pip = pi;
    }
    w(i) = pip * (1.0-pip); // hessian // Unstable as free
    v(i) = y(i) - pi; // gradient
    r(i) = v(i) / w(i); // residual
  }
  // done
  return List::create(Named("beta0")=b0, Named("beta")=beta, Named("r")=r, Named("w")=w);
}


