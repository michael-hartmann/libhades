#ifndef __LIBHADES_ARPACK_H
#define __LIBHADES_ARPACK_H

#include "libhades.h"

void znaupd_(
    int *IDO,         /* reverse communication flag */
    char *BMAT,       /* I: standard, G: generalized eigenproblem */
    int *N,           /* dimension of the eigenproblem */
    char *WHICH,      /* LM, SM, LR, SR, LI, SI */
    int *NEV,         /* number of eigenvalues of OP to be computed. 0 < NEV < N-1. */
    double *TOL,      /* tolerance */
    complex_t *RESID, /* */
    int *NCV,         /* Number of columns of the matrix V */
    complex_t *V,     /* final set of Arnoldi basis vectors */
    int *LDV,         /* Leading dimension of V */
    int *IPARAM,      /* */
    int *IPNTR,       /* */
    complex_t *WORKD, /* */
    complex_t *WORKL, /* */
    int *LWORKL,      /* */
    double *RWORK,    /* */
    int *INFO,        /* */
    int _BMAT,        /* /  The length of the actual BMAT argument */
    int _WHICH        /* The length of the actual WHICH argument */
);

void zneupd_(
    int *RVEC,
    char *HOWMNY,
    int *SELECT, 
    complex_t *D,
    complex_t *Z,
    int *LDZ,
    complex_t *SIGMA,
    complex_t *WORKEV, 
    char *BMAT,
    int *N,
    char *WHICH,
    int *NEV,
    double *TOL, 
    complex_t *RESID,
    int *NCV, 
    complex_t *V,
    int *LDV,
    int *IPARAM,
    int *IPNTR,
    complex_t *WORKD, 
    complex_t *WORKL,
    int *LWORKL,
    double *RWORK,
    int *INFO,
    int _HOWMNY,
    int _BMAT,
    int _WHICH
);

#endif
