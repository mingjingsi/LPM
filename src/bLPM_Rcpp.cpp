// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include "bLPM_aux.hpp"
#include "Update_alpha.hpp"
#include "Update_beta.hpp"
#include "functions.hpp"
#include "cdf_bivnorm.hpp"
using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
RcppExport SEXP bLPM_noX_Rcpp(arma::field<arma::vec>& Pvalue, arma::field<arma::mat>& id, 
                             arma::mat& alpha, arma::mat& pi1_, arma::mat& pair_id,
                             arma::uword maxiter=1e4, double tol=1e-6, const int coreNum=1){
  uword K = Pvalue.n_elem;
  uword npairs = id.n_elem;

  mat fit_alpha = zeros<mat>(2, npairs);
  mat fit_beta0 = zeros<mat>(2, npairs);
  vec fit_rho = zeros<vec>(npairs);

  mat fit_alpha_stage1 = zeros<mat>(2, npairs);
  mat fit_alpha_stage2 = zeros<mat>(2, npairs);
  mat fit_pi1_stage1   = zeros<mat>(2, npairs);
  mat fit_beta0_stage2 = zeros<mat>(2, npairs);

  vec *L_stage1_list[npairs];
  vec *L_stage2_list[npairs];
  vec *L_stage3_list[npairs];

  bLPM_noX bLPM_noXObj(Pvalue, id, K, npairs, pair_id, alpha, pi1_, fit_alpha, fit_alpha_stage1,
                     fit_alpha_stage2, fit_beta0, fit_pi1_stage1, fit_beta0_stage2, fit_rho,
                     maxiter, tol, L_stage1_list, L_stage2_list, L_stage3_list);

  const int n_thread = coreNum;
  std::vector<std::thread> threads(n_thread);
  for(int i_thread = 0; i_thread < n_thread; i_thread++){
    threads[i_thread] = std::thread(&bLPM_noX::fit_by_thread, &bLPM_noXObj, i_thread);
  }

  for(int i = 0; i < n_thread; i++){
    threads[i].join();
  }

  Rcpp::List L_stage1_List(npairs);
  Rcpp::List L_stage2_List(npairs);
  Rcpp::List L_stage3_List(npairs);
  for (int j = 0; j < npairs; j++){
    L_stage1_List[j] = *L_stage1_list[j];
    L_stage2_List[j] = *L_stage2_list[j];
    L_stage3_List[j] = *L_stage3_list[j];
  }

  for(int k = 0; k < npairs; k++){
    delete L_stage1_list[k];
    delete L_stage2_list[k];
    delete L_stage3_list[k];
  }

  Rcpp::List ret;
  ret["alpha"] = bLPM_noXObj.fit_alpha;
  ret["alpha_stage1"] = bLPM_noXObj.fit_alpha_stage1;
  ret["alpha_stage2"] = bLPM_noXObj.fit_alpha_stage2;
  ret["beta0"] = bLPM_noXObj.fit_beta0;
  ret["pi1_stage1"] = bLPM_noXObj.fit_pi1_stage1;
  ret["beta0_stage2"] = bLPM_noXObj.fit_beta0_stage2;
  ret["rho"] = bLPM_noXObj.fit_rho;
  ret["L_stage1_List"] = L_stage1_List;
  ret["L_stage2_List"] = L_stage2_List;
  ret["L_stage3_List"] = L_stage3_List;

  return ret;

}

// [[Rcpp::export]]
RcppExport SEXP bLPM_Rcpp(arma::field<arma::vec>& Pvalue, arma::mat& X, arma::field<arma::mat>& id, 
                        arma::mat& alpha, arma::mat& pi1_, arma::mat& pair_id,
                        arma::uword maxiter=1e4, double tol=1e-6, const int coreNum=1){

  uword K = Pvalue.n_elem;
  uword npairs = id.n_elem;
  uword D = X.n_cols;

  mat fit_alpha = zeros<mat>(2, npairs);
  cube fit_beta = zeros<cube>(2, D, npairs);
  vec fit_rho = zeros<vec>(npairs);

  mat fit_alpha_stage1 = zeros<mat>(2, npairs);
  mat fit_alpha_stage2 = zeros<mat>(2, npairs);
  mat fit_pi1_stage1   = zeros<mat>(2, npairs);
  cube fit_beta_stage2 = zeros<cube>(2, D, npairs);

  vec *L_stage1_list[npairs];
  vec *L_stage2_list[npairs];
  vec *L_stage3_list[npairs];

  bLPM bLPMObj(Pvalue, X, id, K, D, npairs, pair_id, alpha, pi1_, fit_alpha, fit_alpha_stage1,
           fit_alpha_stage2, fit_beta, fit_pi1_stage1, fit_beta_stage2, fit_rho, maxiter, tol,
           L_stage1_list, L_stage2_list, L_stage3_list);

  const int n_thread = coreNum;
  std::vector<std::thread> threads(n_thread);
  for(int i_thread = 0; i_thread < n_thread; i_thread++){
    threads[i_thread] = std::thread(&bLPM::fit_by_thread, &bLPMObj, i_thread);
  }

  for(int i = 0; i < n_thread; i++){
    threads[i].join();
  }

  Rcpp::List L_stage1_List(npairs);
  Rcpp::List L_stage2_List(npairs);
  Rcpp::List L_stage3_List(npairs);
  for (int j = 0; j < npairs; j++){
    L_stage1_List[j] = *L_stage1_list[j];
    L_stage2_List[j] = *L_stage2_list[j];
    L_stage3_List[j] = *L_stage3_list[j];
  }

  for(int k = 0; k < npairs; k++){
    delete L_stage1_list[k];
    delete L_stage2_list[k];
    delete L_stage3_list[k];
  }

  Rcpp::List ret;
  ret["alpha"] = bLPMObj.fit_alpha;
  ret["alpha_stage1"] = bLPMObj.fit_alpha_stage1;
  ret["alpha_stage2"] = bLPMObj.fit_alpha_stage2;
  ret["beta"] = bLPMObj.fit_beta;
  ret["pi1_stage1"] = bLPMObj.fit_pi1_stage1;
  ret["beta_stage2"] = bLPMObj.fit_beta_stage2;
  ret["rho"] = bLPMObj.fit_rho;
  ret["L_stage1_List"] = L_stage1_List;
  ret["L_stage2_List"] = L_stage2_List;
  ret["L_stage3_List"] = L_stage3_List;

  return ret;
}


// [[Rcpp::export]]
RcppExport SEXP single_LPM_noX_Rcpp(arma::vec& Pvalue, double& alpha, double& pi1_, 
                                    arma::uword maxiter, double tol){
  
  uword M = Pvalue.n_rows;
  vec fit_pi1 = zeros<vec>(M);
  vec LL_all = zeros<vec>(maxiter+1);
  uword iter_times = 0;

  Update_alpha_single(Pvalue, M, alpha, pi1_, fit_pi1, iter_times, LL_all, maxiter, tol);

  vec LL = LL_all.subvec(0, iter_times-1);
  
  Rcpp::List ret;
  ret["alpha"] = alpha;
  ret["pi1"] = fit_pi1;
  ret["LL"] = LL;
  
  return ret;
}

// [[Rcpp::export]]
RcppExport SEXP single_LPM_Rcpp(arma::vec& Pvalue, arma::mat& X, double& alpha, double& pi1_,
                                arma::uword maxiter=1e4, double tol=1e-6){

  uword M = Pvalue.n_rows;
  uword D = X.n_cols;

  // Stage 1 (update alpha)
  double alpha_stage1 = alpha;
  vec pi1_stage1 = zeros<vec>(M);
  uword iter_times_stage1 = 0;
  vec L_stage1 = zeros<vec>(maxiter+1);
  Update_alpha_single(Pvalue, M, alpha_stage1, pi1_, pi1_stage1, iter_times_stage1, L_stage1, maxiter, tol);
  vec LL_stage1 = L_stage1.subvec(0, iter_times_stage1-1);
  
  // Stage 2 (update beta)
  double fit_alpha = alpha_stage1;
  vec fit_beta = zeros<vec>(D);
  fit_beta(0) = -qnorm(1 - pi1_);
  vec fit_pi1 = zeros<vec>(M);
  uword iter_times_stage2 = 0;
  vec L_stage2 = zeros<vec>(maxiter+1);
  Update_beta_single(Pvalue, X, M, fit_alpha, fit_beta, fit_pi1, iter_times_stage2, L_stage2, maxiter, tol);
  vec LL_stage2 = L_stage2.subvec(0, iter_times_stage2-1);

  Rcpp::List ret;
  ret["alpha"] = fit_alpha;
  ret["alpha_stage1"] = alpha_stage1;
  ret["beta"] = fit_beta;
  ret["pi1"] = fit_pi1;
  ret["pi1_stage1"] = pi1_stage1;
  ret["L_stage1_List"] = LL_stage1;
  ret["L_stage2_List"] = LL_stage2;

  return ret;
}
