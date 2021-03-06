# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

bLPM_noX_Rcpp <- function(Pvalue, id, alpha, pi1_, pair_id, maxiter = 1e4L, tol = 1e-6, coreNum = 1L) {
    .Call('_LPM_bLPM_noX_Rcpp', PACKAGE = 'LPM', Pvalue, id, alpha, pi1_, pair_id, maxiter, tol, coreNum)
}

bLPM_Rcpp <- function(Pvalue, X, id, alpha, pi1_, pair_id, maxiter = 1e4L, tol = 1e-6, coreNum = 1L) {
    .Call('_LPM_bLPM_Rcpp', PACKAGE = 'LPM', Pvalue, X, id, alpha, pi1_, pair_id, maxiter, tol, coreNum)
}

single_LPM_noX_Rcpp <- function(Pvalue, alpha, pi1_, maxiter, tol) {
    .Call('_LPM_single_LPM_noX_Rcpp', PACKAGE = 'LPM', Pvalue, alpha, pi1_, maxiter, tol)
}

single_LPM_Rcpp <- function(Pvalue, X, alpha, pi1_, maxiter = 1e4L, tol = 1e-6) {
    .Call('_LPM_single_LPM_Rcpp', PACKAGE = 'LPM', Pvalue, X, alpha, pi1_, maxiter, tol)
}

