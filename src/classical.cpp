#include <RcppArmadilloExtensions/sample.h>
#include "helper_functions.h"
using namespace Rcpp;

/*-------------------------------------------------------
# Class to solve for Rzu, M and related quantities given
# observables Rxy, Rxz, Rzy (contained in R) and values
# for K and Rxsu.
#------------------------------------------------------*/
class solveRzu {
public:
  solveRzu(double K, double Rxsu, arma::mat R);
  double Rzu, Su, M, Ruv;
  double Rxy, Rzy, Rxz, A, B1, B2, B, dRxsu, dK;
//private:
};
//Class constructor
solveRzu::solveRzu(double K, double Rxsu, arma::mat R){
  Rxy = R(0, 1);
  Rxz = R(0, 2);
  Rzy = R(1, 2);
  A = Rxsu * Rxz / sqrt(K);
  B1 = (Rxy * Rxz - K * Rzy);
  B2 = sqrt((1 - pow(Rxsu, 2.0)) / (K * (K - pow(Rxy, 2.0))));
  B = B1 * B2;
  Rzu = A - B;
  Su = B1 / (sqrt(K) * Rxz * Rxsu - K * Rzu);
  dRxsu = Rxz/sqrt(K) + Rxsu * B1 /
          sqrt(K * (K - pow(Rxy, 2.0)) * (1 - pow(Rxsu, 2.0)));
  dK = -Rxsu * Rxz / (2 * pow(K, 3.0/2.0)) + B2 *
          (Rzy + B1 / 2.0 * (1 / K + 1/(K - pow(Rxy, 2.0))));
  M = sqrt(1 + pow(dRxsu, 2.0) + pow(dK, 2.0));
  Ruv = (Rxsu * sqrt(K) - Rzu * Rxz) / sqrt(K - pow(Rxz, 2.0));
}




// [[Rcpp::export]]
List classicalSampler(arma::mat Rho_vech, int L, int n_M,
                     double K_L = 0.2, double K_U = 1,
                     double Rxsu_L = -0.9, double Rxsu_U = 0.9,
                     double Rzu_L = -0.9, double Rzu_U = 0.9,
                     double xRsq = 0, bool weight_sampling=true){

  // Rcpp doesn't support arma::cube as function argument
  // so we pass a matrix each of whose columns is the vech
  // of a Rho slice, and then devech them
  int J = Rho_vech.n_cols;
  arma::cube Rho(3, 3, J);
  for(int i = 0; i < J; i++){
    Rho.slice(i) = devech(Rho_vech.col(i), 3);
  }

  int max_iter = n_M * 10;

  Rcpp::IntegerVector indices = Rcpp::seq_len(n_M);
  arma::cube Rzu(L, 1, J, arma::fill::zeros);
  arma::cube Rxsu(L, 1, J, arma::fill::zeros);
  arma::cube K(L, 1, J, arma::fill::zeros);
  arma::cube SuTilde(L, 1, J, arma::fill::zeros);
  arma::cube Ruv(L, 1, J, arma::fill::zeros);
  arma::vec K_lower_bound(J, arma::fill::zeros);
  arma::vec Rxsu_upper_bound(J, arma::fill::zeros);
  arma::vec max_M(J, arma::fill::zeros);
  arma::vec step1_efficiency(J, arma::fill::zeros);

// Loop over draws of observable corr matrix Rho
  for(int j = 0; j < J; j++){
    arma::mat Rhoj = Rho.slice(j);
    arma::vec K_bounds;
    double K_Lower_Inf = (Rhoj(0,1)*Rhoj(0,1)+Rhoj(0,2)*Rhoj(0,2)-2*Rhoj(0,1)*Rhoj(0,2)*Rhoj(1,2))/(1-Rhoj(1,2)*Rhoj(1,2));
    K_bounds << K_L << arma::as_scalar(K_Lower_Inf) << arma::endr;
    double K_Lower = arma::max(K_bounds);
    K_lower_bound(j) = K_Lower;

    // Account for control variables (if present) in bounds for draws of Rxsu
    double RxsuRAW_U = Rxsu_U * sqrt(1 - xRsq / (K_U * (1 - xRsq) + xRsq));
    double Rxsu_Upper = RxsuRAW_U / sqrt(1 - xRsq /
                                           (K_Lower * (1 - xRsq) + xRsq));
    arma::vec Rxsu_bounds;
    Rxsu_bounds << Rxsu_Upper << Rxsu_U << arma::endr;
    Rxsu_Upper = arma::min(Rxsu_bounds);
    Rxsu_upper_bound(j) = Rxsu_Upper;


    // 1st step: reject infeasible points & points outside prior
    int M_count = 0;
    int n_iter = 0;
    arma::vec Rzu_temp(n_M, arma::fill::zeros);
    arma::vec Rxsu_temp(n_M, arma::fill::zeros);
    arma::vec K_temp(n_M, arma::fill::zeros);
    arma::vec M_temp(n_M, arma::fill::zeros);
    arma::vec SuTilde_temp(n_M, arma::fill::zeros);
    arma::vec Ruv_temp(n_M, arma::fill::zeros);

    // Only draw on the identified set if it exists!
    // (User beliefs could be incompatible with Sigma draw)
    if(K_U > K_Lower){

      while((M_count < n_M) && (n_iter < max_iter)){

        double K_sim = R::runif(K_Lower, K_U);
        double Rxsu_sim = R::runif(Rxsu_L, Rxsu_U);
        solveRzu soln(K_sim, Rxsu_sim, Rhoj);
        bool Ruv_constraint = (pow(soln.Ruv, 2.0) +
                               pow(soln.Rzu, 2.0)) < 1.0;
        bool Rzu_constraint = (soln.Rzu < Rzu_U) &&
                              (soln.Rzu > Rzu_L);

        if(Rzu_constraint && Ruv_constraint){
          Rzu_temp(M_count) = soln.Rzu;
          M_temp(M_count) = soln.M;
          Rxsu_temp(M_count) = Rxsu_sim;
          K_temp(M_count) = K_sim;
          SuTilde_temp(M_count) = soln.Su;
          Ruv_temp(M_count) = soln.Ruv;
          M_count++;
        }// end if
        n_iter++;
      }// end while

      // 2nd step: weighted resampling
      if(n_iter < max_iter){

        if(weight_sampling){

          //resampling probabilities
          double max_Mj  = arma::max(M_temp);
          arma::vec pM = M_temp / max_Mj;
          arma::uvec rand_indices = as<arma::uvec>(
            RcppArmadillo::sample(indices, L, true,
                                  NumericVector(pM.begin(), pM.end())));
          rand_indices = rand_indices - 1; //arma uses zero-indexing!
          Rzu.slice(j) = Rzu_temp(rand_indices);
          Rxsu.slice(j) = Rxsu_temp(rand_indices);
          K.slice(j) = K_temp(rand_indices);
          SuTilde.slice(j) = SuTilde_temp(rand_indices);
          Ruv.slice(j) = Ruv_temp(rand_indices);
          max_M(j) = max_Mj;
          step1_efficiency(j) = double(n_M) / double(n_iter);

        } else {

          // Sampling without weights (Take first L draws)
          Rcpp::IntegerVector rand_indices_list = Rcpp::seq_len(L);
          arma::uvec rand_indices = as<arma::uvec>(rand_indices_list)-1;

          double max_Mj  = arma::max(M_temp);
          Rzu.slice(j) = Rzu_temp(rand_indices);
          Rxsu.slice(j) = Rxsu_temp(rand_indices);
          K.slice(j) = K_temp(rand_indices);
          SuTilde.slice(j) = SuTilde_temp(rand_indices);
          Ruv.slice(j) = Ruv_temp(rand_indices);
          max_M(j) = max_Mj;
          step1_efficiency(j) = double(n_M) / double(n_iter);

        }


      }// end if
    }// end if

  }// end for

  return List::create(Named("Rzu") = Rzu,
                      Named("Rxsu") = Rxsu,
                      Named("K") = K,
                      Named("SuTilde") = SuTilde,
                      Named("Ruv") = Ruv,
                      Named("Klower") = K_lower_bound,
                      Named("Rxsu_Upper") = Rxsu_upper_bound,
                      Named("step1eff") = step1_efficiency,
                      Named("maxM") = max_M);
}
