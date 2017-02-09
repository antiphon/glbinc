#include <Rcpp.h>
#include "general.h"
using namespace Rcpp;


// Cross product of the jth column of x with y
double crossprod(NumericMatrix &x, NumericVector& y, int l) {
  double val = 0;
  for (int i=0; i< y.size(); i++) val += x(i,l) * y(i);
  return(val);
}

void gd_binomial(NumericMatrix& beta,
                 NumericMatrix &X,
                 NumericVector& r,
                 NumericVector& eta,
                 int g,
                 IntegerVector& G1,
                 int l,
                 double l1, double l2,
                 NumericVector& df,
                 NumericVector& betaold){
  int G, i, j;
  int n = X.nrow();
  double nd = (double) n;
  double shift, si;
  double len, v = 0.25;

  G = G1(g+1) - G1(g);// group size
  NumericVector z(G);
  double z_norm;

  for(j=G1(g); j < G1(g+1); j++) z(j-G1(g)) = crossprod(X, r, j)/nd + betaold(j);
  z_norm = norm(z);
  len = S(v * z_norm, l1) / ( v * (1 + l2));

  if(len != 0 | betaold(G1(g)) != 0) {
    for(j = G1(g); j < G1(g+1); j++) {
      beta(j,l) = len * z(j-G1(g)) / z_norm;
      shift = beta(j,l) - betaold(j);
      for(i = 0; i < n; i++) {
        si = shift * X(i,j);
        r(i) -= si;
        eta(i) += si;
      }
    }
  }
  if(len > 0) df(l) += G * len / z_norm;
}

// [[Rcpp::export]]

List glbin_lcd_cpp(NumericMatrix X,
                   NumericVector y,
                   NumericVector offset,
                   int G0,
                   IntegerVector G1,
                   NumericVector group_weight,
                   NumericVector lambda,
                   bool add_intercept,
                   double alpha,
                   int verb,
                   double eps,
                   int dfmax,
                   int maxiter,
                   int AIC_stop
) {

  int n = X.nrow();
  int p = X.ncol();

  double nd = (double) n;

  int i, j, g, it, violations, converged;

  int L = lambda.size();
  int J = G1.size() - 1;
  // Outcome without intercept
  NumericMatrix beta(p, L);
  NumericVector b0(L);
  NumericVector df(L);
  NumericVector AIC(L);
  NumericVector AIC_check(AIC_stop + 1);

  if(verb>3) {
    Rprintf("Got:\n");
    if(add_intercept) Rprintf("* intercept\n");
    else Rprintf("* no intercept\n");
    Rprintf("* unpenalised terms: %i\n", G0);
    Rprintf("* penalised groups: ");
    for(g = 0; g < J; g++) Rprintf("% i", G1(g+1)-G1(g));
    Rprintf("\n* eps: %f\n* alpha: %f\n* maxiter: %i\n* dfmax %i\n* AIC stop: %i\n", eps, alpha, maxiter, dfmax, AIC_stop);
  }

  //
  NumericVector betaold(p);
  IntegerVector e(J,0);
  double ybar = mean(y);
  double a0;
  b0(0) = log(ybar/(1.0-ybar));
  a0 = log(ybar/(1.0-ybar));
  NumericVector lik(L);
  NumericVector eta(n);
  NumericVector r(n);
  for(i=0; i < n; i++){
    eta(i) = a0 + offset(i);
  }
  //a0 = a0 + offset(0);
  double aic_old = 9999999;
  char aic_sign[] = "-";

  double pi, pip, shift, si, l1, l2;
  int ndf=0;
  // Path
  for(int l = 0; l < L; l++) {
    Rcpp::checkUserInterrupt();

    if(l>0){
      a0 = b0(l-1);
      for(j = 0; j < p; j++) betaold(j) = beta(j, l-1);
    }


    // check if df above dfmax
    if(ndf > dfmax) {
      // stop
      if(verb) Rprintf("dfmax reached.\n");
      break;
    }

    // check if AIC increasing
    if(AIC_stop > 0) if(l > AIC_stop){
      if(verb>100) Rprintf("AIC check:");
      for(j = 0; j <= AIC_stop; j++) {
        AIC_check(j) = AIC(l - AIC_stop + j - 1);
        if(verb>100)Rprintf(" %f", AIC_check(j));
      }
      if(verb>100) Rprintf("\n");
      // now check each consecutive is increasing
      i = 0;
      for(j=0; j < AIC_stop; j++){
        if((AIC_check(j+1) - AIC_check(j)) > 0) i++;
      }
      if(i==AIC_stop) {
        if(verb) Rprintf("AIC increasing.\n");
        break;
      }
    }
    // update linear predictor
    // for(i=0; i < n; i++){
    //   eta(i) = a0 + offset(i);
    //   for(j = 0; j < p; j++){
    //     eta(i) += X(i,j) * betaold(j);
    //   }
    // }

    // check dfmax

    it = 0;
    double unpenalised_shifts;
    while(it < maxiter) {
      while(it < maxiter) {
        converged = 0;
        it++;
        if(verb) Rprintf("        \r[%i/%i it:%i]", l+1, L, it);
        // quadratic approximation of the log likelihood
        lik(l) = 0;
        for(i=0; i < n; i++) {
          if(eta(i) > 10) {
            pi = 1;
            pip = 1-1e-5;
          }else if(eta(i) < -10){
            pi = 0;
            pip = 1e-5;
          } else{
            pi = 1.0/(1.0 + exp(-eta(i)));
            pip = pi;
          }
          r(i) = (y(i)-pi) / 0.25;
          lik(l) += y(i) * log(pip) + (1-y(i)) * log(1-pip);
        }

        df(l) = 0;
        unpenalised_shifts = 0;
        // update intercept if needed
        if(add_intercept){
          shift = mean(r);
          b0(l) = shift + a0;
          for(i=0; i < n; i++) {
            r[i] -= shift;
            eta[i] += shift;
          }
          df(l)++;
          unpenalised_shifts += std::abs(shift);
        }

        // update unpenalised groups
        for(j=0; j < G0; j++) {
          shift = crossprod(X, r, j)/nd;
          unpenalised_shifts += std::abs(shift);
          beta(j,l) = shift + betaold(j);
          for(i=0; i < n; i++) {
            si = shift * X(i,j);
            r(i) -= si;
            eta(i) += si;
          }
          df(l)++;
        }

        // update penalised terms
        for(g = 0; g < J; g++) {
          if(e(g)!=0) {
            l1 = lambda(l) * group_weight(g) * alpha;// lasso
            l2 = lambda(l) * group_weight(g) * (1-alpha);// ridge
            gd_binomial(beta, X, r, eta, g, G1, l, l1, l2, df, betaold);
          }
        }

        // check convergence
        if(unpenalised_shifts < eps){
          if(checkConvergence(beta, betaold, l, eps)){
          converged++;
          break;
          }
        }
        a0 = b0(l);
        for(j = 0; j < p; j++) betaold(j) = beta(j,l);
      } // inner while loop

      // scan for newly activated groups
      violations=0;
      for(g = 0; g < J; g++) {
        if(e(g)==0) {
          l1 = lambda(l) * group_weight(g) * alpha;// lasso
          l2 = lambda(l) * group_weight(g) * (1-alpha);// ridge
          gd_binomial(beta, X, r, eta, g, G1, l, l1, l2, df, betaold);
          if(beta(G1(g),l) != 0) { // check if group has activated
            e(g) = 1;
            violations++;
          }
        }
      }
      if(violations == 0) break;
      a0 = b0(l);
      for(j = 0; j < p; j++) betaold(j) = beta(j,l);
    } // middle loop

    // check df
    ndf = G0 + 1;
    for(g=0; g < J; g++) {
      if(betaold(G1(g)) != 0)  ndf += G1(g+1)-G1(g);
    }
    // gather AIC
    AIC(l) = -2*lik(l) + 2*df(l);
    aic_sign[0] = (aic_old > AIC(l)) ? '-' : '+';
    aic_old = AIC(l);
    if(verb) Rprintf(" [lambda %f][non-zero %i rel_df %f][aic %s]\n", lambda[l], ndf, df(l), aic_sign);

  } // outer loop over the path

  // done
  return List::create(Named("beta0")=b0, Named("beta")=beta, Named("df")=df, Named("lik")=lik, Named("AIC")=AIC);
}



