#ifndef bLPM_aux_hpp
#define bLPM_aux_hpp

#include <stdio.h>
#include <RcppArmadillo.h>
#include <thread>
#include <mutex>

using namespace std;
using namespace arma;

class bLPM_noX{
public:
  uword current_idx=0;

  field<vec> Pvalue;
  field<mat> id;
  uword K;
  uword npairs;
  mat pair_id;
  mat alpha;
  mat pi1_;
  mat fit_alpha;
  mat fit_alpha_stage1;
  mat fit_alpha_stage2;
  mat fit_beta0;
  mat fit_pi1_stage1;
  mat fit_beta0_stage2;
  vec fit_rho;
  uword maxiter;
  double tol;
  vec **L_stage1_list;
  vec **L_stage2_list;
  vec **L_stage3_list;

  bLPM_noX(const field<vec>& Pvalue, const field<mat>& id, const uword K, const uword npairs, 
          const mat& pair_id, const mat& alpha,
         const mat& pi1_, mat& fit_alpha, mat& fit_alpha_stage1, mat& fit_alpha_stage2,
         mat& fit_beta0, mat& fit_pi1_stage1, mat& fit_beta0_stage2, vec& fit_rho,
         uword maxiter, double tol, vec *L_stage1_list[], vec *L_stage2_list[],
         vec *L_stage3_list[]){

    this -> Pvalue = Pvalue;
    this -> id = id;
    this -> K = K;
    this -> npairs = npairs;
    this -> pair_id = pair_id;
    this -> alpha = alpha;
    this -> pi1_ = pi1_;
    this -> fit_alpha = fit_alpha;
    this -> fit_alpha_stage1 = fit_alpha_stage1;
    this -> fit_alpha_stage2 = fit_alpha_stage2;
    this -> fit_beta0 = fit_beta0;
    this -> fit_pi1_stage1 = fit_pi1_stage1;
    this -> fit_beta0_stage2 = fit_beta0_stage2;
    this -> fit_rho = fit_rho;
    this -> maxiter = maxiter;
    this -> tol = tol;
    this -> L_stage1_list = L_stage1_list;
    this -> L_stage2_list = L_stage2_list;
    this -> L_stage3_list = L_stage3_list;

  }

  void loop_by_thread(int i);
  void fit_by_thread(int thread_id);
  int  next();
};

class bLPM{
public:
  uword current_idx=0;

  field<vec> Pvalue;
  mat X;
  field<mat> id;
  uword K;
  uword D;
  uword npairs;
  mat pair_id;
  mat alpha;
  mat pi1_;
  mat fit_alpha;
  mat fit_alpha_stage1;
  mat fit_alpha_stage2;
  cube fit_beta;
  mat fit_pi1_stage1;
  cube fit_beta_stage2;
  vec fit_rho;
  uword maxiter;
  double tol;
  vec **L_stage1_list;
  vec **L_stage2_list;
  vec **L_stage3_list;

  bLPM(const field<vec>& Pvalue, const mat& X, const field<mat>& id, const uword K, const uword D,
     const uword npairs, const mat& pair_id, const mat& alpha, const mat& pi1_, mat& fit_alpha,
     mat& fit_alpha_stage1, mat& fit_alpha_stage2, cube& fit_beta, mat& fit_pi1_stage1,
     cube& fit_beta_stage2, vec& fit_rho, uword maxiter, double tol, vec *L_stage1_list[],
     vec *L_stage2_list[], vec *L_stage3_list[]){

    this -> Pvalue = Pvalue;
    this -> X = X;
    this -> id = id;
    this -> K = K;
    this -> D = D;
    this -> npairs = npairs;
    this -> pair_id = pair_id;
    this -> alpha = alpha;
    this -> pi1_ = pi1_;
    this -> fit_alpha = fit_alpha;
    this -> fit_alpha_stage1 = fit_alpha_stage1;
    this -> fit_alpha_stage2 = fit_alpha_stage2;
    this -> fit_beta = fit_beta;
    this -> fit_pi1_stage1 = fit_pi1_stage1;
    this -> fit_beta_stage2 = fit_beta_stage2;
    this -> fit_rho = fit_rho;
    this -> maxiter = maxiter;
    this -> tol = tol;
    this -> L_stage1_list = L_stage1_list;
    this -> L_stage2_list = L_stage2_list;
    this -> L_stage3_list = L_stage3_list;
  }

  void loop_by_thread(int i);
  void fit_by_thread(int thread_id);
  int  next();
};

#endif

