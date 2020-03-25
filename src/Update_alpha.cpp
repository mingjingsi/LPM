// EM algorithm for Stage 1

#include "Update_alpha.hpp"
#include "functions.hpp"

void Update_alpha(const mat& Pvalue, const uword M, vec& alpha, vec& pi1_, mat& pi1,
                  uword& iter_times, vec& LL, const uword maxiter, const double tol){

  double L = 0;
  double L_old = L;

  mat PvalueT = Pvalue.t();

  for (uword iter = 0; iter < maxiter; iter++){
    // E step
    mat comp_pos = pi1_.t() % alpha.t() % pow_mat_vec(Pvalue, alpha - 1).each_row();
    mat comp_pos1 = zeros<mat>(M, 2);
    comp_pos1.col(0) = comp_pos.col(0) + 1 - pi1_(0);
    comp_pos1.col(1) = comp_pos.col(1) + 1 - pi1_(1);
    pi1 = comp_pos/comp_pos1;

    // compute incomplete log likelihood
    L = accu(log(comp_pos1));
    LL(iter) = L;

    // M step
    mat pi1T = pi1.t();
    vec colsum_pi1 = sum(pi1T, 1);
    alpha = -colsum_pi1/sum(pi1T%log(PvalueT), 1);
    pi1_ = colsum_pi1/M;

    // check convergence
    if (iter != 0){
      if (L < L_old){
        printf("L is not increasing. \n");
        break;
      }
      if ((L - L_old)/abs(L) < tol){
        iter_times = iter + 1;
        break;
      }
      else{
        iter_times = maxiter;
      }
    }
    L_old = L;
  }
}


void Update_alpha_single(const vec& Pvalue, const uword M, double& alpha, double& pi1_, 
                         vec& pi1, uword& iter_times, vec& LL, const uword maxiter, const double tol){
  
  double L = 0;
  double L_old = L;

  for (uword iter = 0; iter < maxiter; iter++){
    // E step
    vec comp_pos = (pi1_*alpha) * pow(Pvalue,(alpha-1));
    pi1 = comp_pos/(comp_pos+1-pi1_);
    
    // compute incomplete log likelihood
    L = accu(log(comp_pos+1-pi1_));
    LL(iter) = L;
    
    // M step
    alpha = -sum(pi1)/dot(pi1,log(Pvalue));
    pi1_ = sum(pi1)/M;
    
    // check convergence
    if (iter != 0){
      if (L < L_old){
        printf("L is not increasing. \n");
        break;
      }
      if ((L - L_old)/abs(L) < tol){
        iter_times = iter + 1;
        break;
      }
      else{
        iter_times = maxiter;
      }
    }
    L_old = L;
  }
}

