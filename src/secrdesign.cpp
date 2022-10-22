/*
   External procedures for secrdesign package 2022-10-21

*/

#include "secrdesign.h"

/*==============================================================================*/

arma::mat hazmat (const arma::vec par, const arma::mat &d, int detectfn) {
    
    arma::mat H = d;
    
    int kk = d.n_rows;    /* number of traps */
    int mm = d.n_cols;    /* number of points on mask */
    int k,m;

    //-------------------------------------------------------------------------
    
    if (detectfn == 0 || detectfn == 14) {       // HHN
        H = par(0) * arma::exp(-arma::square(d) / 2 / par(1) / par(1));
    }
    else if (detectfn == 1 || detectfn == 15) {  // HHR
        H = par(0) * (1 - arma::exp(- arma::pow(d /par(1), - par(2))));
    } 
    else if (detectfn == 2 || detectfn == 16) {  // HEX
        H = par(0) * arma::exp(-d / par(1));
    } 
    else if (detectfn == 6 || detectfn == 17) {  // HAN
        H = par(0) * exp(-arma::square(d-par(2)) / 2 / par(1) / par(1));
    } 
    else if (detectfn == 8 || detectfn == 18) {  // HCG
        for (k=0; k<kk; k++) {
            for (m=0; m<mm; m++) {
                boost::math::gamma_distribution<> gam(par(2), par(1)/par(2));
                H(k,m) = par(0) * boost::math::cdf(complement(gam,d(k,m)));
            }
        }
    } 
    else if (detectfn == 19) {                   // HVP
        H = par(0) * arma::exp(- arma::pow(d / par(1), par(2)));
    }
    else Rcpp::stop ("detectfn not implemented");
    
    return (H);
}
/*==============================================================================*/

// [[Rcpp::export]]
Rcpp::List Lambdacpp (int type, const arma::vec par, const arma::mat d, int detectfn)
{
    // traps x mask hazard matrix
    arma::mat h = hazmat(par, d, detectfn);
    
    arma::rowvec sumhk = arma::sum(h, 0); 
    Rcpp::NumericVector outsumhk = Rcpp::NumericVector(sumhk.begin(), sumhk.end());

    Rcpp::NumericVector outsumpk (1);
    Rcpp::NumericVector outsumq2 (1);
    outsumpk[0] = NA_REAL;
    outsumq2[0] = NA_REAL;
    
    arma::rowvec sumhk2 = arma::sum(arma::square(h), 0); 
    arma::rowvec sumq2 = sumhk2 / arma::square(sumhk);	// sumhk>0
    outsumq2 = Rcpp::NumericVector(sumq2.begin(), sumq2.end());
    
    // multi
    if (type == 0) {
        arma::rowvec sumpk = 1 - arma::exp(- sumhk);
        outsumpk = Rcpp::NumericVector(sumpk.begin(), sumpk.end());
    }
    // proximity
    else if (type == 1) {
        arma::rowvec sumpk = sum(1 - arma::exp(- h), 0);
        outsumpk = Rcpp::NumericVector(sumpk.begin(), sumpk.end());
    }
    // count: sumpk not used
    
    return Rcpp::List::create(
        Rcpp::Named("sumhk") = outsumhk,
        Rcpp::Named("sumpk") = outsumpk,
        Rcpp::Named("sumq2") = outsumq2
    );
}
/*==============================================================================*/

// [[Rcpp::export]]
Rcpp::List Qpmcpp (
        const arma::vec par, 
        const arma::rowvec D, 
        const arma::mat d, 
        int detectfn,
        int noccasions)
{
    double Qp, Qpm;
    double G = arma::accu(D);  // assume already includes cellsize

    // traps x mask matrix hazard per occasion
    arma::mat h = hazmat(par, d, detectfn);

    // mask probability detection
    arma::rowvec pd = 1 - arma::exp(- arma::sum(h,0) * noccasions);
    arma::rowvec p0 = 1-pd;
    pd = D % pd;
    Qp = arma::accu (pd) / G;

    arma::mat pij = 1 - arma::exp(- h * noccasions);
    pij = pij/(1-pij);
    arma::rowvec p1 = p0 % arma::sum(pij, 0);   // sum over traps
    arma::rowvec p2 = D % (1 - p0 - p1);
    Qpm = arma::accu (p2) / G;

    return (Rcpp::List::create(
            Rcpp::Named("Qp") = Qp,
            Rcpp::Named("Qpm") = Qpm));
}
/*==============================================================================*/
