#include "Rcpp.h"
using namespace Rcpp;
// 2025-06-09 slower than R!

// [[Rcpp::export]]
NumericMatrix simOUcpp (
        NumericVector xy, 
        double beta, 
        double sigma, 
        int noccasions, 
        NumericVector start) {
    NumericMatrix out (noccasions, 2);
    NumericVector mean (2);
    int i,j;
    double sd = sigma * sqrt(1 - exp(-2*beta));
    out(1,_) = start;
    for (i=1; i<noccasions; i++) {
        mean = (1 - exp(-beta)) * xy + exp(-beta) * out(i - 1,_);
        for (j=0;j<2;j++) 
            out(i,j) = R::rnorm(mean(j), sd);
    }
    return out;
}
