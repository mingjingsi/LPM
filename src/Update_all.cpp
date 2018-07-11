// PX-EM algorithm for Stage 3

#include "Update_all.hpp"
#include "cdf_bivnorm.hpp"
#include "functions.hpp"

void Update_all_noX(const mat& Pvalue, const uword M, vec& alpha, vec& beta, double& rho,
                    vec& pi11, vec& pi10, vec& pi01, vec& pi00, uword& iter_times,
                    vec& LL, const uword maxiter, const double tol){

  double L = 0;
  double L_old = L;

for (uword iter = 0; iter < maxiter; iter++){
    // E step

    // E_eta
    vec Phi = pnorm(beta);
    double Phi11 = bivnor(-beta(0), -beta(1), rho);
    double Phi10 = -Phi11 + Phi(0);
    double Phi01 = -Phi11 + Phi(1);
    double Phi00 = 1 + Phi11 - Phi(0) - Phi(1);

    mat pow_p_alpha = pow_mat_vec(Pvalue, alpha - 1);
    mat alpha_pow = alpha.t()%pow_p_alpha.each_row();

    vec alpha_pow11 = alpha_pow.col(0)%alpha_pow.col(1);

    vec comp_pi11 = Phi11*alpha_pow11;
    vec comp_pi10 = Phi10*alpha_pow.col(0);
    vec comp_pi01 = Phi01*alpha_pow.col(1);
    double comp_pi00 = Phi00;
    vec comp_pi   = comp_pi11 + comp_pi10 + comp_pi01 + comp_pi00;

    pi11 = comp_pi11/comp_pi;
    pi10 = comp_pi10/comp_pi;
    pi01 = comp_pi01/comp_pi;
    pi00 = comp_pi00/comp_pi;

    vec E_eta1 = pi11 + pi10;
    vec E_eta2 = pi11 + pi01;


    // E_Z
    double c = 1/sqrt(1 - rho*rho);

    vec phi = dnorm(beta);

    double red = phi(0)*gsl_cdf_ugaussian_P((beta(1) - rho*beta(0))*c);
    double green = phi(1)*gsl_cdf_ugaussian_P((beta(0) - rho*beta(1))*c);

    double comp1 = red + rho*green;
    double comp2 = green + rho*red;
    double comp3 = -beta(0)*red - rho*rho*beta(1)*green +
        rho/c*phi(1)*gsl_ran_gaussian_pdf((beta(0) - rho*beta(1))*c, 1);
    double comp4 = -beta(1)*green - rho*rho*beta(0)*red +
        rho/c*phi(0)*gsl_ran_gaussian_pdf((beta(1) - rho*beta(0))*c, 1);
    double comp5 = -rho*beta(0)*red - rho*beta(1)*green +
        phi(0)/c*gsl_ran_gaussian_pdf((beta(1) - rho*beta(0))*c, 1);

    vec pi_Phi = (alpha_pow11 -alpha_pow.col(0) - alpha_pow.col(1) + 1)/comp_pi;
    vec pi_Phi10_00_phi1 = (alpha_pow.col(0) - 1)*phi(0)/comp_pi;
    vec pi_Phi01_00_phi2 = (alpha_pow.col(1) - 1)*phi(1)/comp_pi;

    vec E_Z1 = beta(0) + pi_Phi*comp1 + pi_Phi10_00_phi1 + rho*pi_Phi01_00_phi2;
    vec E_Z2 = beta(1) + pi_Phi*comp2 + rho*pi_Phi10_00_phi1 + pi_Phi01_00_phi2;

    vec E_Z1_2 = beta(0)*(2*E_Z1 - beta(0)) + 1 + pi_Phi*comp3 -
      beta(0)*pi_Phi10_00_phi1 - rho*rho*beta(1)*pi_Phi01_00_phi2;
    vec E_Z2_2 = beta(1)*(2*E_Z2 - beta(1)) + 1 + pi_Phi*comp4 -
      rho*rho*beta(0)*pi_Phi10_00_phi1 - beta(1)*pi_Phi01_00_phi2;
    vec E_Z1_Z2 = -beta(0)*beta(1) + beta(0)*E_Z2 + beta(1)*E_Z1 + rho + pi_Phi*comp5 -
      rho*beta(0)*pi_Phi10_00_phi1 - rho*beta(1)*pi_Phi01_00_phi2;

    // compute incomplete log likelihood
    L = accu(log(comp_pi));
    LL(iter) = L;

    // M step
    alpha(0) = -sum(E_eta1)/sum(E_eta1%log(Pvalue.col(0)));
    alpha(1) = -sum(E_eta2)/sum(E_eta2%log(Pvalue.col(1)));

    double mean_E_Z1 = mean(E_Z1);
    double mean_E_Z2 = mean(E_Z2);
    mat Sigma = zeros<mat>(2, 2);
    Sigma(0, 0) = mean(E_Z1_2) - mean_E_Z1*mean_E_Z1;
    Sigma(0, 1) = mean(E_Z1_Z2) - mean_E_Z1*mean_E_Z2;
    Sigma(1, 0) = Sigma(0, 1);
    Sigma(1, 1) = mean(E_Z2_2) - mean_E_Z2*mean_E_Z2;

    beta(0) = mean_E_Z1/sqrt(Sigma(0, 0));
    beta(1) = mean_E_Z2/sqrt(Sigma(1, 1));
    rho = Sigma(0, 1)/sqrt(Sigma(0, 0))/sqrt(Sigma(1, 1));

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


void Update_all(const mat& Pvalue, const mat& X, const uword M, vec& alpha, mat& beta,
                double& rho, vec& pi11, vec& pi10, vec& pi01, vec& pi00, uword& iter_times,
                vec& LL, const uword maxiter, const double tol){

  double L = 0;
  double L_old = L;

  mat XT = X.t();
  mat XTX = XT*X;
  mat XTX_1 = inv(XTX);

  for (uword iter = 0; iter < maxiter; iter++){
    // E step

    // E_eta
    mat Xbeta = X * beta.t();
    mat Phi = pnorm_mat(Xbeta);
    vec Phi11 = zeros<vec>(M);
    for (uword i = 0; i < M; i++)
      Phi11(i) = bivnor(-Xbeta(i, 0), -Xbeta(i, 1), rho);

    vec Phi10 = -Phi11 + Phi.col(0);
    vec Phi01 = -Phi11 + Phi.col(1);
    vec Phi00 = 1 + Phi11 - Phi.col(0) - Phi.col(1);

    mat pow_p_alpha = pow_mat_vec(Pvalue, alpha - 1);
    mat alpha_pow = alpha.t()%pow_p_alpha.each_row();

    vec alpha_pow11 = alpha_pow.col(0)%alpha_pow.col(1);

    vec comp_pi11 = Phi11%alpha_pow11;
    vec comp_pi10 = Phi10%alpha_pow.col(0);
    vec comp_pi01 = Phi01%alpha_pow.col(1);
    vec comp_pi00 = Phi00;
    vec comp_pi   = comp_pi11 + comp_pi10 + comp_pi01 + comp_pi00;

    pi11 = comp_pi11/comp_pi;
    pi10 = comp_pi10/comp_pi;
    pi01 = comp_pi01/comp_pi;
    pi00 = comp_pi00/comp_pi;

    vec E_eta1 = pi11 + pi10;
    vec E_eta2 = pi11 + pi01;

    // E_Z
    double c = 1/sqrt(1 - rho*rho);

    mat phi = dnorm_mat(Xbeta);

    vec red = phi.col(0)%pnorm((Xbeta.col(1) - rho*Xbeta.col(0))*c);
    vec green = phi.col(1)%pnorm((Xbeta.col(0) - rho*Xbeta.col(1))*c);

    vec Xbeta1_red = Xbeta.col(0)%red;
    vec Xbeta2_green = Xbeta.col(1)%green;
    vec phi1_dnorm = phi.col(0)%dnorm((Xbeta.col(1) - rho*Xbeta.col(0))*c);

    vec comp1 = red + rho*green;
    vec comp2 = green + rho*red;
    vec comp3 = -Xbeta1_red - rho*rho*Xbeta2_green +
      rho/c*phi.col(1)%dnorm((Xbeta.col(0) - rho*Xbeta.col(1))*c);
    vec comp4 = -Xbeta2_green - rho*rho*Xbeta1_red + rho/c*phi1_dnorm;
    vec comp5 = -rho*Xbeta1_red - rho*Xbeta2_green + phi1_dnorm/c;

    vec pi_Phi = (alpha_pow11 -alpha_pow.col(0) - alpha_pow.col(1) + 1)/comp_pi;
    vec pi_Phi10_00_phi1 = (alpha_pow.col(0) - 1)%phi.col(0)/comp_pi;
    vec pi_Phi01_00_phi2 = (alpha_pow.col(1) - 1)%phi.col(1)/comp_pi;

    vec E_Z1 = Xbeta.col(0) + pi_Phi%comp1 + pi_Phi10_00_phi1 + rho*pi_Phi01_00_phi2;
    vec E_Z2 = Xbeta.col(1) + pi_Phi%comp2 + rho*pi_Phi10_00_phi1 + pi_Phi01_00_phi2;

    mat E_Z = zeros<mat>(M, 2);
    E_Z.col(0) = E_Z1;
    E_Z.col(1) = E_Z2;

    vec pi_Phi10_00_phi1_Xbeta1 = pi_Phi10_00_phi1%Xbeta.col(0);
    vec pi_Phi01_00_phi2_Xbeta2 = pi_Phi01_00_phi2%Xbeta.col(1);

    vec E_Z1_2 = Xbeta.col(0)%(2*E_Z1 - Xbeta.col(0)) + 1 + pi_Phi%comp3 -
      pi_Phi10_00_phi1_Xbeta1 - rho*rho*pi_Phi01_00_phi2_Xbeta2;
    vec E_Z2_2 = Xbeta.col(1)%(2*E_Z2 - Xbeta.col(1)) + 1 + pi_Phi%comp4 -
      rho*rho*pi_Phi10_00_phi1_Xbeta1 - pi_Phi01_00_phi2_Xbeta2;
    vec E_Z1_Z2 = -Xbeta.col(0)%Xbeta.col(1) + Xbeta.col(0)%E_Z2 + Xbeta.col(1)%E_Z1 +
      rho + pi_Phi%comp5 - rho*pi_Phi10_00_phi1_Xbeta1 - rho*pi_Phi01_00_phi2_Xbeta2;

    // compute incomplete log likelihood
    L = accu(log(comp_pi));
    LL(iter) = L;

    // M step
    alpha(0) = -sum(E_eta1)/sum(E_eta1%log(Pvalue.col(0)));
    alpha(1) = -sum(E_eta2)/sum(E_eta2%log(Pvalue.col(1)));

    mat XT_E_Z = XT*E_Z;
    mat gamma = XT_E_Z.t()*XTX_1;
    mat sum_EZZ = zeros<mat>(2, 2);
    sum_EZZ(0, 0) = sum(E_Z1_2);
    sum_EZZ(1, 0) = sum(E_Z1_Z2);
    sum_EZZ(0, 1) = sum_EZZ(1, 0);
    sum_EZZ(1, 1) = sum(E_Z2_2);

    mat Sigma = (sum_EZZ - gamma*XT_E_Z)/M;
    beta = gamma.each_col()/sqrt(Sigma.diag());
    rho = Sigma(0, 1)/sqrt(Sigma(0, 0))/sqrt(Sigma(1, 1));

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
