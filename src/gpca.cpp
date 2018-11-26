#define ARMA_64BIT_WORD 1
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// ugmd = matrix(nrow = n,ncol = k)
//   vgmd = matrix(nrow = p, ncol = k)
//   dgmd = rep(0,k)
//   propv = rep(0,k)
//   Xhat = X
//   
//   u = rnorm(n)
//   v = rnorm(p)
//   
//   qrnorm = sum(diag(t(X) %*% Q %*% X %*% R))
//   
//   cumv = rep(0,k)
//   
//   thr = 1e-6
// for(i in 1:k){
//   err=1
//   while(err > thr){
//     oldu = u
//     oldv = v
//     uhat = Xhat %*% R %*% v
//     u = uhat/as.double(sqrt(t(uhat)%*% Q %*% uhat))
//     vhat = t(Xhat) %*% Q %*% u
//     v = vhat/as.double(sqrt(t(vhat) %*% R %*% vhat))
//     err = t(oldu - u) %*% (oldu - u) + t(oldv -v ) %*% (oldv - v)
//   }
//   dgmd[i] = t(u) %*% Q %*% X %*% R %*% v
//     ugmd[,i] = u
//     vgmd[,i] = v
//     Xhat = Xhat - dgmd[i] *  u %*% t(v)
//     propv[i] = dgmd[i]^2/as.double(qrnorm)
//     cumv[i] = sum(propv[1:i])
// }


// large dimension should be in ROWS (if columns is X, then pass (t(X), R,))

//[[Rcpp::export]]
List gmd_deflation_cpp(const arma::mat &X, arma::sp_mat Q, arma::sp_mat R, int k, double thr=1e-5) {
  Rcout << "begin " << std::endl;
  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat ugmd(n, k, fill::zeros);
  arma::mat vgmd(p, k, fill::zeros);
  arma::vec dgmd(k, fill::zeros);
  arma::vec propv(k, fill::zeros);
  arma::vec cumv(k, fill::zeros);
  Rcout << "assign Xhat " << std::endl;
  arma::mat Xhat = X;
  
  arma::vec u = randn(n);
  arma::vec v = randn(p);
  
  Rcout << "qrnorm " << std::endl;
  
  double qrnorm = 1;
  
  if (p > n) {
    qrnorm = trace(X * R * X.t() * Q);
  } else {
    qrnorm = trace(X.t() * Q * X * R);
  }
  
  Rcout << "begin loop " << std::endl;
  for (int i=0; i<k; i++) {
    double err=1;
    while (err > thr) {
      arma::vec oldu = vec(u);
      arma::vec oldv = vec(v);
      
      arma::vec uhat = Xhat * R * v;
      u = uhat/sqrt(uhat.t() * Q * uhat).at(0,0);
      arma::vec vhat = Xhat.t() * Q * u;
      v = vhat/sqrt(vhat.t() * R * vhat).at(0,0);
      
      err = sum(((oldu - u).t() * (oldu - u)) + ((oldv -v ).t() * (oldv - v)));
      //Rcout << "error: " << err << std::endl;
      
    }
    
    dgmd(i) = (u.t() * Q * X * R * v).eval().at(0,0);
    ugmd.col(i) = u;
    vgmd.col(i) = v;
    Xhat = Xhat - dgmd(i) *  u * v.t();
    propv(i) = pow(dgmd(i),2)/qrnorm;
    cumv(i) = sum(propv);
    //dgmd(i) = tmp(0,0);
    //dgmd(i) = (u.t() * Q * X * R * v)
  }
  
  return List::create(Named("d")=dgmd, Named("v")=vgmd, Named("u")=ugmd, Named("cumv")=cumv, Named("propv")=propv);
  
}

//[[Rcpp::export]]
List gmdLA_cpp(const arma::mat &X, arma::sp_mat Q, arma::mat R, int k) {
  int n = X.n_rows;
  //int p = X.n_cols;
  arma::vec r_eigval;
  arma::mat r_eigvec;
  arma::eig_sym(r_eigval, r_eigvec, R);
  
  arma::uvec keep = arma::reverse(arma::find(r_eigval > 0 && arma::abs(r_eigval) > 1e-7));
  r_eigvec = r_eigvec.cols(keep);
  r_eigval = r_eigval(keep);
  
  arma::mat Rtilde = r_eigvec * arma::diagmat(sqrt(r_eigval)) * r_eigvec.t();
  
  arma::vec inv_values = 1 / sqrt(r_eigval);
  arma::mat Rtilde_inv = r_eigvec * arma::diagmat(inv_values) * r_eigvec.t();
  
  //return List::create(r_eigvec, r_eigval, keep, Rtilde);

  arma::mat inmat = X.t() * Q * X;
  arma::mat RtinRt = Rtilde * inmat * Rtilde;
  arma::mat XR = X * R;
  // 
  
  arma::mat RnR = R * inmat * R;
  arma::vec eigval;
  arma::mat eigvec;
   
  arma::eig_sym(eigval, eigvec, RtinRt);
  
  keep = arma::reverse(arma::find(arma::abs(eigval) > 1e-7));
  k = std::min<int>(keep.n_elem, k);

  arma::mat vgmd = Rtilde_inv * eigvec.cols(keep.subvec(0,k-1));
  arma::vec dgmd = sqrt(eigval(keep.subvec(0,k-1)));
  arma::mat ugmd(n,k);
  
  //propv <- dgmd ^ 2 / sum(diag(as.matrix(inmat %*% R)))
  double normalizing_number = 1;
  
  for (int i=0; i<k; i++) {
   // normalizing.number = sqrt(vgmd[, i] %*% RnR %*% vgmd[, i])
    normalizing_number = sqrt(vgmd.col(i).t() * RnR * vgmd.col(i)).at(0,0);
    ugmd.col(i) = XR * vgmd.col(i) / normalizing_number;
  }
  
  // cumv = rep(0, k)
  // propv <- dgmd ^ 2 / sum(diag(as.matrix(inmat %*% R)))
  // normalizing.number <- 1
  // 
  return List::create(Named("d")=dgmd, Named("u")=ugmd, Named("v")=vgmd);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
library(Matrix)
#Q <- sparseMatrix(i=c(2,5,7,9), j=c(5,7,12,15), x=c(2,2,2,2), dims=c(15,15))
#R <- sparseMatrix(i=c(2,5,7,9), j=c(5,7,12,15), x=c(3,2,2,2),dims=c(15,15))

N=100
M=500000
X <- matrix(rnorm(N*M),N,M)
Q <- neighborweights::spatial_smoother(as.matrix(1:N), nnk=12) + diag(x=1, nrow(X), nrow(X))
R <- neighborweights::spatial_smoother(as.matrix(1:M), nnk=12) + Diagonal(x=1, n=ncol(X))
gmd_deflation_cpp(X,Q,R,5)
gmd_deflation_cpp(t(X),R,Q,5)

gmdLA_cpp(X,Q,R,10)

#gmdLA_cpp(X, Q, as.matrix(R), 1,N,N) 

*/
