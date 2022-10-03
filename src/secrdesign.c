/*
   External procedures for secrdesign package

   can compile with gcc 4.6.3 :
   gcc -Ic:/R/R-3.3.0/include -c secrdesign.c -Wall -pedantic -std=gnu99

*/

/* 2017-03-26 created */
/* 2017-08-02 L2 output from sumpkC */
/* 2017-10-31 detectfn HAN, HCG bug fixed */
/* 2020-11-02 protect sumpkC against divide by zero */
/* 2022-04-11 LambdaK for capped detectors */

#include "secrdesign.h"
#include <time.h>

/*==============================================================================*/

void R_CheckUserInterrupt(void);

/*==============================================================================*/

double d2 (
        int k,
        int m,
        double A1[],
                 double A2[],
                          int A1rows,
                          int A2rows)
    /*
    return squared distance between two points given by row k in A1
    and row m in A2, where A1 and A2 have respectively A1rows and A2rows
    */
{
    return(
        (A1[k] - A2[m]) * (A1[k] - A2[m]) +
            (A1[k + A1rows] - A2[m + A2rows]) * (A1[k + A1rows] - A2[m + A2rows])
    );
}

void LambdaC (
        double *par,       /* lambda0, sigma, z */
    int    *kk,        /* number of traps */
    int    *mm,        /* number of points on mask */
    double *traps,     /* x,y locations of traps (first x, then y) */
    double *mask,      /* x,y points on mask (first x, then y) */
    int    *fn,        /* detectfn code 14 = halfnormal */
    double *L,         /* return value vector of length mm */
    int    *resultcode /* 0 for successful completion */
)
{
    int k,m;
    double d2km;
    *resultcode = 1;                   /* generic failure */
    
    for (m=0; m<*mm; m++) {
        L[m] = 0;
        for (k = 0; k < *kk; k++) {
            d2km = d2(k, m, traps, mask, *kk, *mm);
            if (*fn == 14)
                L[m] += exp(-d2km / 2 / par[1] / par[1]);
            else if (*fn == 15)
                L[m] += 1 - exp(- pow(sqrt(d2km) / par[1] , - par[2]));
            else if (*fn == 16)
                L[m] += exp(-sqrt(d2km) / par[1]);
            else if (*fn == 17)
                L[m] += exp(- sqrt(d2km-par[2]) * sqrt(d2km-par[2]) / 2 / par[1] / par[1]);
            else if (*fn == 18)
                L[m] += pgamma(sqrt(d2km), par[2], par[1]/par[2], 0,0);
            else error("only detectfn 14:18");
        }
        L[m] *= par[0];
    }
    *resultcode = 0;                   /* successful completion */
}

void LambdaK (
    double *par,       /* lambda0, sigma, z */
    int    *kk,        /* number of traps */
    int    *mm,        /* number of points on mask */
    double *traps,     /* x,y locations of traps (first x, then y) */
    double *mask,      /* x,y points on mask (first x, then y) */
    int    *fn,        /* detectfn code 14 = halfnormal */
    double *LK,        /* return value vector of length kk */
    int    *resultcode /* 0 for successful completion */
)
{
    int k,m;
    double d2km;
    *resultcode = 1;                   /* generic failure */
    
    for (k = 0; k < *kk; k++) {
        LK[k] = 0;
        for (m=0; m<*mm; m++) {
            d2km = d2(k, m, traps, mask, *kk, *mm);
            if (*fn == 14)
                LK[k] += exp(-d2km / 2 / par[1] / par[1]);
            else if (*fn == 15)
                LK[k] += 1 - exp(- pow(sqrt(d2km) / par[1] , - par[2]));
            else if (*fn == 16)
                LK[k] += exp(-sqrt(d2km) / par[1]);
            else if (*fn == 17)
                LK[k] += exp(- sqrt(d2km-par[2]) * sqrt(d2km-par[2]) / 2 / par[1] / par[1]);
            else if (*fn == 18)
                LK[k] += pgamma(sqrt(d2km), par[2], par[1]/par[2], 0,0);
            else error("only detectfn 14:18");
        }
        LK[k] *= par[0];
        
        
    }
    *resultcode = 0;                   /* successful completion */
}

void Encapped (
        double *par,       /* lambda0, sigma, z */
    int    *kk,        /* number of traps */
    int    *mm,        /* number of points on mask */
    double *traps,     /* x,y locations of traps (first x, then y) */
    double *mask,      /* x,y points on mask (first x, then y) */
    int    *fn,        /* detectfn code 0 = halfnormal */
    int    *nocc,
    double *LK,         /* LK from LambdaK; includes D factor*/
    double *D,         /* density per cell */
    double *En,         /* result */
    int    *resultcode /* 0 for successful completion */
)
{
    int k,m;
    double d2km;
    double lam, prod, ps;
    *resultcode = 1;                   /* generic failure */
    
    *En = 0;
    for (m=0; m<*mm; m++) {
        prod = 1.0;
        for (k = 0; k < *kk; k++) {
            lam = par[0];
            d2km = d2(k, m, traps, mask, *kk, *mm);
            if (*fn == 14)
                lam *= exp(-d2km / 2 / par[1] / par[1]);
            else if (*fn == 15)
                lam *= 1 - exp(- pow(sqrt(d2km) / par[1] , - par[2]));
            else if (*fn == 16)
                lam *= exp(-sqrt(d2km) / par[1]);
            else if (*fn == 17)
                lam *= exp(- sqrt(d2km-par[2]) * sqrt(d2km-par[2]) / 2 / par[1] / par[1]);
            else if (*fn == 18)
                lam *= pgamma(sqrt(d2km), par[2], par[1]/par[2], 0,0);
            else error("only detectfn 14:18");
            
            ps = (1 - exp(-LK[k])) * lam/LK[k];
            prod *= 1 - ps;   // not detected on one occasion
        }
        *En += (1 - pow(prod, *nocc)) * *D;
        
    }
    *resultcode = 0;                   /* successful completion */
}


void sumpkC (
        int    *type,      /* 0 multi 1 proximity*/
    double *par,       /* lambda0, sigma, z */
    int    *kk,        /* number of traps */
    int    *mm,        /* number of points on mask */
    double *traps,     /* x,y locations of traps (first x, then y) */
    double *mask,      /* x,y points on mask (first x, then y) */
    int    *fn,        /* detectfn code 14 = halfnormal */
    double *L,         /* return value vector of length mm */
    double *L2,        /* return value vector of length mm */
    int    *resultcode /* 0 for successful completion */
)
{
    int k,m;
    double hk, pk, ps, sumhk, sumhk2, d2km, dkm;
    
    *resultcode = 1;                   /* generic failure */
    if (*type > 2) error("unrecognised type in sumpkC");
    for (m=0; m<*mm; m++) {
        L[m] = 0;
        L2[m] = 0;
        sumhk = 0;
        sumhk2 = 0;
        for (k = 0; k < *kk; k++) {
            d2km = d2(k, m, traps, mask, *kk, *mm);
            dkm = sqrt(d2km);
            if (*fn == 14)
                hk = par[0] * exp(-d2km / 2 / par[1] / par[1]);
            else if (*fn == 15)
                hk = par[0] * (1 - exp(- pow(sqrt(d2km) / par[1] , - par[2])));
            else if (*fn == 16)
                hk = par[0] * exp(-dkm / par[1]);
            else if (*fn == 17) 
                hk = par[0] * exp(- (dkm-par[2]) * (dkm-par[2])  / 2 / par[1] / par[1]);	    
            else if (*fn == 18)
                hk = par[0] * pgamma(dkm, par[2], par[1]/par[2], 0,0);
            else error("only detectfn 14:18");
            
            if (*type == 1) {
                pk = 1 - exp(- hk);
                L[m] += pk;
            }
            sumhk += hk;
            sumhk2 += hk*hk;
        }
        if (*type == 0) {
            ps =  1 - exp(-sumhk);
            L[m] = ps;
        }
        if (sumhk>0) {   // 2020-11-02
            L2[m] = sumhk2 / sumhk / sumhk;	
        }
    }
    *resultcode = 0;                   /* successful completion */
}
/*==============================================================================*/

// UNDER DEVELOPMENT 2017-05-06

void repeatr (
        int    *type,      /* 0 multi 1 proximity 2 count */
double *par,       /* lambda0, sigma, z */
    int    *kk,        /* number of traps */
    int    *mm,        /* number of points on mask */
    double *traps,     /* x,y locations of traps (first x, then y) */
    double *mask,      /* x,y points on mask (first x, then y) */
    int    *fn,        /* detectfn code 14 = halfnormal */
    double *L,         /* return value vector of length mm */
    int    *resultcode /* 0 for successful completion */
)
{
    int k, m;
    double sumhk = 0;
    double p1 = 0;
    double * h = NULL;
    if (*type > 2) error("unrecognised type in repeatr");
    h = (double *)  S_alloc(*kk, sizeof (double));
    *resultcode = 1;                   /* generic failure */
    if (*fn != 14) error("only hazard halfnormal");
    for (m=0; m<*mm; m++) {
	L[m] = 0;
	sumhk = 0;
	for (k = 0; k < *kk; k++) {
	    h[k] = par[0] * exp(-d2(k, m, traps, mask, *kk, *mm) /2 / par[1] / par[1]);
	    if (*type == 0) 
		sumhk += h[k];
	}
	for (k = 0; k < *kk; k++) {
	    // p1 is Pr(ct at k given caught somewhere)
	    if (*type == 0) 
		p1 = h[k] / (1 - exp(-sumhk)) / (1 - exp(-sumhk));
	    else if (*type == 1 || *type == 2) 
		p1 = (1 - exp(- h[k] )) / (1 - exp(-sumhk));
	    L[m] += p1 * p1; 
	}
	L[m] = 1 - L[m];
    }
    *resultcode = 0;                   /* successful completion */
}
/*==============================================================================*/