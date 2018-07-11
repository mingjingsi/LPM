#include "functions.hpp"
using namespace arma;

vec pnorm(vec x){
  uword M = x.n_elem;
  vec Phi = zeros<vec>(M);
  for (uword i = 0; i < M; i++)
    Phi(i) = gsl_cdf_ugaussian_P(x(i));

  return Phi;
}

mat pnorm_mat(mat x){
  uword M = x.n_rows;
  uword N = x.n_cols;
  mat Phi = zeros<mat>(M, N);
  for (uword i = 0; i < M; i++){
    for (uword j = 0; j < N; j++)
      Phi(i, j) = gsl_cdf_ugaussian_P(x(i, j));
  }

  return Phi;
}

vec dnorm(vec x){
  uword M = x.n_elem;
  vec phi = zeros<vec>(M);
  for (uword i = 0; i < M; i++)
    phi(i) = gsl_ran_gaussian_pdf(x(i), 1);

  return phi;
}

mat dnorm_mat(mat x){
  uword M = x.n_rows;
  uword N = x.n_cols;
  mat phi = zeros<mat>(M, N);
  for (uword i = 0; i < M; i++){
    for (uword j = 0; j < N; j++)
      phi(i, j) = gsl_ran_gaussian_pdf(x(i, j), 1);
  }

  return phi;
}

vec qnorm(vec x){
  uword M = x.n_elem;
  vec q = zeros<vec>(M);
  for (uword i = 0; i < M; i++)
    q(i) = gsl_cdf_ugaussian_Pinv(x(i));

  return q;
}

mat pow_mat_vec(mat x, vec y){
  uword M = x.n_rows;
  mat pow_x_y = zeros<mat>(M, 2);
  pow_x_y.col(0) = pow(x.col(0), y(0));
  pow_x_y.col(1) = pow(x.col(1), y(1));

  return(pow_x_y);
}
