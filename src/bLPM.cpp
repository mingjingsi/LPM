#include "bLPM_aux.hpp"
#include "Update_alpha.hpp"
#include "Update_beta.hpp"
#include "Update_all.hpp"
#include "functions.hpp"

using namespace std;
using namespace arma;

void bLPM_noX::loop_by_thread(int i){

  uword current_id1 = pair_id(0, i) - 1;
  uword current_id2 = pair_id(1, i) - 1;
  uword current_M = id(i).n_rows;
  mat current_Pvalue = zeros<mat>(current_M, 2);
  for (uword j = 0; j < current_M; j++){
    current_Pvalue(j, 0) = Pvalue(current_id1)(id(i)(j, 0) - 1);
    current_Pvalue(j, 1) = Pvalue(current_id2)(id(i)(j, 1) - 1);
  }

  printf("current pair %d\n", i+1);

  // Stage 1 (update alpha)
  vec alpha_stage1 = alpha.col(i);
  vec pi1_1 = pi1_.col(i);
  mat pi1_stage1 = zeros<mat>(current_M, 2);
  uword iter_times_stage1 = 0;
  vec L_stage1 = zeros<vec>(maxiter+1);
  Update_alpha(current_Pvalue, current_M, alpha_stage1, pi1_1, pi1_stage1, iter_times_stage1, L_stage1,
               maxiter, tol);
  vec LL_stage1 = L_stage1.subvec(0, iter_times_stage1-1);

  // Stage 2 (update beta)
  vec alpha_stage2 = alpha_stage1;
  vec beta0_stage2 = -qnorm(1 - pi1_1);
  mat pi1_stage2 = zeros<mat>(current_M, 2);
  uword iter_times_stage2 = 0;
  vec L_stage2 = zeros<vec>(maxiter+1);
  Update_beta_noX(current_Pvalue, current_M, alpha_stage2, beta0_stage2, pi1_stage2, iter_times_stage2,
                  L_stage2, maxiter, tol);
  vec LL_stage2 = L_stage2.subvec(0, iter_times_stage2-1);

  // Stage 3 (update all)
  vec alpha_stage3 = alpha_stage2;
  vec beta0_stage3 = beta0_stage2;
  double rho = 0;
  vec pi11 = zeros<vec>(current_M);
  vec pi10 = zeros<vec>(current_M);
  vec pi01 = zeros<vec>(current_M);
  vec pi00 = zeros<vec>(current_M);
  uword iter_times_stage3 = 0;
  vec L_stage3 = zeros<vec>(maxiter+1);
  Update_all_noX(current_Pvalue, current_M, alpha_stage3, beta0_stage3, rho, pi11, pi10, pi01, pi00,
                 iter_times_stage3, L_stage3, maxiter, tol);
  vec LL_stage3 = L_stage3.subvec(0, iter_times_stage3-1);

  fit_alpha.col(i) = alpha_stage3;
  fit_beta0.col(i) = beta0_stage3;
  fit_rho(i) = rho;

  fit_alpha_stage1.col(i) = alpha_stage1;
  fit_alpha_stage2.col(i) = alpha_stage2;
  fit_pi1_stage1.col(i)   = pi1_1;
  fit_beta0_stage2.col(i) = beta0_stage2;

  L_stage1_list[i] = new vec(LL_stage1.begin(),LL_stage1.n_elem, true);
  L_stage2_list[i] = new vec(LL_stage2.begin(),LL_stage2.n_elem, true);
  L_stage3_list[i] = new vec(LL_stage3.begin(),LL_stage3.n_elem, true);

}

std::mutex mtx;
int bLPM_noX::next(){
  std::lock_guard<std::mutex> lock(mtx);
  if(current_idx >= npairs){
    return -1;
  }
  current_idx++;
  return current_idx-1;
}

void bLPM_noX::fit_by_thread(int thread_id){
  while(true){
    int idx = next();

    if(idx == -1){
      break;
    }

    loop_by_thread(idx);
  }
}

void bLPM::loop_by_thread(int i){

  uword current_id1 = pair_id(0, i) - 1;
  uword current_id2 = pair_id(1, i) - 1;
  uword current_M = id(i).n_rows;
  mat current_Pvalue = zeros<mat>(current_M, 2);
  mat current_X = zeros<mat>(current_M, D);
  for (uword j = 0; j < current_M; j++){
    current_Pvalue(j, 0) = Pvalue(current_id1)(id(i)(j, 0) - 1);
    current_Pvalue(j, 1) = Pvalue(current_id2)(id(i)(j, 1) - 1);
    current_X.row(j) = X.row(id(i)(j, 2) - 1);
  }

  printf("current pair %d\n", i+1);

  // Stage 1 (update alpha)
  vec alpha_stage1 = alpha.col(i);
  vec pi1_1 = pi1_.col(i);
  mat pi1_stage1 = zeros<mat>(current_M, 2);
  uword iter_times_stage1 = 0;
  vec L_stage1 = zeros<vec>(maxiter+1);
  Update_alpha(current_Pvalue, current_M, alpha_stage1, pi1_1, pi1_stage1, iter_times_stage1, L_stage1,
               maxiter, tol);
  vec LL_stage1 = L_stage1.subvec(0, iter_times_stage1-1);

  // Stage 2 (update beta)
  vec alpha_stage2 = alpha_stage1;
  mat beta_stage2 = zeros<mat>(2, D);
  beta_stage2.col(0) = -(qnorm(1 - pi1_1));
  mat pi1_stage2 = zeros<mat>(current_M, 2);
  uword iter_times_stage2 = 0;
  vec L_stage2 = zeros<vec>(maxiter+1);
  Update_beta(current_Pvalue, current_X, current_M, alpha_stage2, beta_stage2, pi1_stage2, iter_times_stage2,
              L_stage2, maxiter, tol);
  vec LL_stage2 = L_stage2.subvec(0, iter_times_stage2-1);

  // Stage 3 (update all)
  vec alpha_stage3 = alpha_stage2;
  mat beta_stage3 = beta_stage2;
  double rho = 0;
  vec pi11 = zeros<vec>(current_M);
  vec pi10 = zeros<vec>(current_M);
  vec pi01 = zeros<vec>(current_M);
  vec pi00 = zeros<vec>(current_M);
  uword iter_times_stage3 = 0;
  vec L_stage3 = zeros<vec>(maxiter+1);
  Update_all(current_Pvalue, current_X, current_M, alpha_stage3, beta_stage3, rho, pi11, pi10, pi01, pi00,
             iter_times_stage3, L_stage3, maxiter, tol);
  vec LL_stage3 = L_stage3.subvec(0, iter_times_stage3-1);

  fit_alpha.col(i) = alpha_stage3;
  fit_beta.slice(i) = beta_stage3;
  fit_rho(i) = rho;

  fit_alpha_stage1.col(i) = alpha_stage1;
  fit_alpha_stage2.col(i) = alpha_stage2;
  fit_pi1_stage1.col(i)   = pi1_1;
  fit_beta_stage2.slice(i) = beta_stage2;

  L_stage1_list[i] = new vec(LL_stage1.begin(),LL_stage1.n_elem, true);
  L_stage2_list[i] = new vec(LL_stage2.begin(),LL_stage2.n_elem, true);
  L_stage3_list[i] = new vec(LL_stage3.begin(),LL_stage3.n_elem, true);

}

// std::mutex mtx;
int bLPM::next(){
  std::lock_guard<std::mutex> lock(mtx);
  if(current_idx >= npairs){
    return -1;
  }
  current_idx++;
  return current_idx-1;

}

void bLPM::fit_by_thread(int thread_id){
  while(true){
    int idx = next();

    if(idx == -1){
      break;
    }
    loop_by_thread(idx);
  }
}
