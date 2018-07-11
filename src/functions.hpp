#ifndef functions_hpp
#define functions_hpp
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

using namespace arma;

vec pnorm(vec x);

mat pnorm_mat(mat x);

vec dnorm(vec x);

mat dnorm_mat(mat x);

vec qnorm(vec x);

mat pow_mat_vec(mat x, vec y);


#endif


