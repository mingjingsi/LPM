#ifndef Update_beta_hpp
#define Update_beta_hpp
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
using namespace arma;

void Update_beta_noX(const mat& Pvalue, const uword M, vec& alpha, vec& beta, mat& pi1,
                     uword& iter_times, vec& LL, const uword maxiter, const double tol);

void Update_beta(const mat& Pvalue, const mat& X, const uword M, vec& alpha, mat& beta,
                 mat& pi1, uword& iter_times, vec& LL, const uword maxiter, const double tol);

void Update_beta_single(const vec& Pvalue, const mat& X, const uword M, double& alpha, vec& beta,
                        vec& pi1, uword& iter_times, vec& LL, const uword maxiter, const double tol);

#endif
