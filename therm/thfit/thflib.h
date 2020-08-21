#ifndef _THFLIB_H_
#define _THFLIB_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>

#define THF_RU 8.3144621
#define THF_CAL2J 4.184

#define THF_NLST 7
#define THF_NHOP 5
#define THF_NEXT 231
#define THF_NCOF 14
#define THF_RWRK 2587
#define THF_IWRK 32
#define THF_MAX1L 256

#if defined(_WIN32) && !defined(__CYGWIN32__) && !defined(__MINGW32__)
#define MINPACK_LMDIF LMDIF
#define LAPACK_DGESDD DGESDD
#else
#define MINPACK_LMDIF lmdif_
#define LAPACK_DGESDD dgesdd_
#endif

#ifdef __cplusplus
extern "C" {
#endif
void MINPACK_LMDIF(void (*fcn)(int *m, int *n, double *x, double *fvec, int *iflag),
                   int *m, int *n, double *x, double *fvec, double *ftol, double *xtol,
                   double *gtol, int *maxfev, double *epsfcn, double *diag,
                   int *mode, double *factor, int *nprint, int *info, int *nfev,
                   double *fjac, int *ldfjac, int *ipvt, double *qtf,
                   double *wa1, double *wa2, double *wa3, double *wa4);
void LAPACK_DGESDD(char *jobz, int *m, int *n, double **a,
                   int *lda, double *s, double **u, int *ldu, double **vt, int *ldvt,
                   double *work, int *lwork, int *iwork, int *info);
#ifdef __cplusplus
}
#endif

typedef struct {
  double *temp;
  double *cp;
  double *weight;
  int ndata;
  int natoms;
  int rotn;
} cpdat_t;

double cpVib(double tdT);
double cpElec(double a, double edT1, double edT2);
double cpHOE1(double T, double *p);
double cpHOE2(double T, double *p);
double cpHOEpoly(double T, double *p, int natoms, int rotn);
double cpHOE(double T, double *p, int natoms, int rotn);

void cpFcn(int *m, int *n, double *x, double *fvec, int *iflag);
int cpExtrap(double *params, double *temp, double *cp, double *weight, int ndata, int natoms,
             int rotn, double tol, int maxfev, int *iwa, double *wa);

double nasa7Cp(double *coeffs, double T);
double nasa7CpD1(double *coeffs, double T);
double nasa7CpD2(double *coeffs, double T);
double nasa7H(double *coeffs, double T);
double nasa7S(double *coeffs, double T);

int cpFit1(double *temp, double *cp, double *weight, double tempbase, double cpbase, int ndata,
           double *coeffs, double *rwork, int *iwork);
int cpFit(double *temp, double *cp, int ndata, int itmin, int itmax, int itoverlap, double woverlap,
          double h298, double s298, double *coeffs, double *rwork, int *iwork);


# endif
