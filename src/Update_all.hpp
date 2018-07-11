#ifndef Update_all_hpp
#define Update_all_hpp
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
using namespace arma;

void Update_all_noX(const mat& Pvalue, const uword M, vec& alpha, vec& beta, double& rho,
                    vec& pi11, vec& pi10, vec& pi01, vec& pi00, uword& iter_times,
                    vec& LL, const uword maxiter, const double tol);

void Update_all(const mat& Pvalue, const mat& X, const uword M, vec& alpha, mat& beta,
                double& rho, vec& pi11, vec& pi10, vec& pi01, vec& pi00, uword& iter_times,
                vec& LL, const uword maxiter, const double tol);

#endif
