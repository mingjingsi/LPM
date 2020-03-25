#ifndef Update_alpha_hpp
#define Update_alpha_hpp
#include <RcppArmadillo.h>
using namespace arma;

void Update_alpha(const mat& Pvalue, const uword M, vec& alpha, vec& pi1_, mat& pi1,
                  uword& iter_times, vec& LL, const uword maxiter, const double tol);

void Update_alpha_single(const vec& Pvalue, const uword M, double& alpha, double& pi1_, 
                         vec& pi1, uword& iter_times, vec& LL, const uword maxiter, const double tol);

#endif
