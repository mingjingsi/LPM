// PX-EM algorithm for Stage 2

#include "Update_beta.hpp"
#include "functions.hpp"

void Update_beta_noX(const mat& Pvalue, const uword M, vec& alpha, vec& beta, mat& pi1,
                     uword& iter_times, vec& LL, const uword maxiter, const double tol){

  double L = 0;
  double L_old = L;
  mat PvalueT = Pvalue.t();

  for (uword iter = 0; iter < maxiter; iter++){
    // E step
    vec Phi = pnorm(beta);
    vec phi = dnorm(beta);

    mat pow_p_alpha = pow_mat_vec(Pvalue, alpha - 1);
    mat comp_pi1 = Phi.t() % alpha.t() % pow_p_alpha.each_row();

    vec comp_pi0 = 1 - Phi;
    mat comp_L = comp_pi1.each_row() + comp_pi0.t();
    pi1 = comp_pi1/comp_L;

    mat v = phi.t() % alpha.t() % pow_p_alpha.each_row()/comp_L -
      phi.t() /comp_L.each_row();
    mat E_Z = beta.t() + v.each_row();
    mat E_Z2 = beta.t() % E_Z.each_row() + 1;

    // compute incomplete log likelihood
    L = accu(log(comp_L));
    LL(iter) = L;

    // M step
    mat pi1T = pi1.t();
    vec gamma = mean(E_Z.t(), 1);
    vec sigma2 = mean(E_Z2.t(), 1) - square(gamma);
    beta = gamma/sqrt(sigma2);
    alpha = -sum(pi1T, 1)/sum(pi1T%log(PvalueT), 1);

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

void Update_beta(const mat& Pvalue, const mat& X, const uword M, vec& alpha, mat& beta,
                 mat& pi1, uword& iter_times, vec& LL, const uword maxiter, const double tol){

  double L = 0;
  double L_old = L;

  mat XTX_1 = inv(X.t()*X);
  mat PvalueT = Pvalue.t();

  for (uword iter = 0; iter < maxiter; iter++){
    // E step
    mat Xbeta = X * beta.t();
    mat Phi = pnorm_mat(Xbeta);
    mat phi = dnorm_mat(Xbeta);

    mat pow_p_alpha = pow_mat_vec(Pvalue, alpha - 1);
    mat Phi_pow = Phi % pow_p_alpha;
    mat comp_pi1 = alpha.t() % Phi_pow.each_row();

    mat comp_pi0 = 1 - Phi;
    mat comp_L = comp_pi1 + comp_pi0;
    pi1 = comp_pi1/comp_L;

    mat phi_pow = phi % pow_p_alpha;
    mat v = alpha.t() % phi_pow.each_row()/comp_L - phi/comp_L;
    mat E_Z = Xbeta + v;
    mat E_Z2 = Xbeta % E_Z + 1;

    // compute incomplete log likelihood
    L = accu(log(comp_L));
    LL(iter) = L;

    // M step
    mat E_ZT = E_Z.t();
    mat pi1T = pi1.t();
    mat gamma = E_ZT*X*XTX_1;
    mat Xgamma = X*gamma.t();
    mat XgammaT = Xgamma.t();
    vec sigma2 = mean(E_Z2.t() - 2*E_ZT%XgammaT + XgammaT%XgammaT, 1);
    beta = gamma.each_col()/sqrt(sigma2);
    alpha = -sum(pi1T, 1)/sum(pi1T%log(PvalueT), 1);

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


void Update_beta_single(const vec& Pvalue, const mat& X, const uword M, double& alpha, vec& beta,
                 vec& pi1, uword& iter_times, vec& LL, const uword maxiter, const double tol){
  
  double L = 0;
  double L_old = L;
  
  mat XTX_1 = inv(X.t()*X);

  for (uword iter = 0; iter < maxiter; iter++){
    // E step
    vec Xbeta = X * beta;
    vec Phi = pnorm(Xbeta);
    vec phi = dnorm(Xbeta);
    
    vec pow_p_alpha = pow(Pvalue, alpha - 1);
    vec comp_pi1 = alpha * Phi % pow_p_alpha;
    
    vec comp_pi0 = 1 - Phi;
    vec comp_L = comp_pi1 + comp_pi0;
    pi1 = comp_pi1/comp_L;
    
    vec v = alpha * phi % pow_p_alpha/comp_L - phi /comp_L;
    vec E_Z = Xbeta + v;
    vec E_Z2 = Xbeta % E_Z + 1;
    
    // compute incomplete log likelihood
    L = accu(log(comp_L));
    LL(iter) = L;
    
    // M step
    vec gamma = (E_Z.t()*X*XTX_1).t();
    vec Xgamma = X*gamma;
    double sigma2 = mean(E_Z2 - 2*E_Z%Xgamma + Xgamma%Xgamma);
    beta = gamma/sqrt(sigma2);
    alpha = -sum(pi1)/sum(pi1%log(Pvalue));
    
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
