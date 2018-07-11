// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include "bLPM_aux.hpp"
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

