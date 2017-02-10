#include <RcppArmadillo.h>
#include <stdexcept>
using namespace Rcpp;

//' Simulate draws from the inverse Wishart distribution
//' @param n An integer, the number of draws.
//' @param v An integer, the degrees of freedom of the distribution.
//' @param S A numeric matrix, the scale matrix of the distribution.
//' @return A numeric array of matrices, each of which is one simulation draw.
//' @details Employs the Bartlett Decomposition (Smith & Hocking 1972).
//' Output exactly matches that of riwish from the MCMCpack package if the same
//' random seed is used.
// [[Rcpp::export]]
arma::cube rinvwish(int n, int v, arma::mat S){
  RNGScope scope;
  int p = S.n_rows;
  arma::mat L = arma::chol(arma::inv_sympd(S), "lower");
  arma::cube sims(p, p, n, arma::fill::zeros);
  for(int j = 0; j < n; j++){
  arma::mat A(p,p, arma::fill::zeros);
  for(int i = 0; i < p; i++){
    int df = v - (i + 1) + 1; //zero-indexing
    A(i,i) = sqrt(R::rchisq(df));
  }
  for(int row = 1; row < p; row++){
    for(int col = 0; col < row; col++){
      A(row, col) = R::rnorm(0,1);
    }
  }
  arma::mat LA_inv = arma::inv(arma::trimatl(arma::trimatl(L) *
    arma::trimatl(A)));
  sims.slice(j) = LA_inv.t() * LA_inv;
  }
  return(sims);
}
