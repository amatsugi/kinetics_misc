#include "thflib.h"

cpdat_t cpdat;

double cpVib(double tdT) {
  double exmt = exp(-tdT);
  return tdT*tdT * exmt / (exmt-1.) / (exmt-1.);
}
double cpElec(double a, double edT1, double edT2) {
  // preliminary implemention
  double tmp = 1. + a*exp(-edT1) + (1.-a)*exp(-edT1);
  double tmp1 = (a*edT1*edT1*exp(-edT1) + (1.-a)*edT2*edT2*exp(-edT2)) / tmp;
  double tmp2 = (a*edT1*exp(-edT1) + (1.-a)*edT2*exp(-edT2)) / tmp;
  return tmp1 - tmp2*tmp2;
}

double cpHOE1(double T, double *p) {
  // electronic contrib.
  double a=p[0], te1=p[1], te2=p[2];
  return (2.5 + cpElec(a, te1/T, te2/T)) * THF_RU;
}

double cpHOE2(double T, double *p) {
  // vib. contrib. (single HO) only
  // TODO: include electronic contrib. if necessary
  double t1=p[0];
  return (3.5 + cpVib(t1/T)) * THF_RU;
}

double cpHOEpoly(double T, double *p, int natoms, int rotn) {
  // vib. contrib. (3 HO model of Ritter1991)
  // TODO: include electronic contrib. if necessary
  double b1=p[0], b2=p[1], t1=p[2], t2=p[3], t3=p[4];
  double cvib;
  int nvib = 3 * natoms - 6 - rotn;
  cvib = b1*cpVib(t1/T) + b2*cpVib(t2/T) + (nvib - (b1+b2)) * cpVib(t3/T);
  return (4. + rotn/2. + cvib) * THF_RU;
}

double cpHOE(double T, double *p, int natoms, int rotn) {
  if (natoms == 1) { return cpHOE1(T, p); }
  else if (natoms == 2) { return cpHOE2(T, p); }
  else { return cpHOEpoly(T, p, natoms, rotn); }
}

void cpFcn(int *m, int *n, double *x, double *fvec, int *iflag) {
  int i;
  for (i=0; i<cpdat.ndata; i++) {
    fvec[i] = cpdat.weight[i] * (cpdat.cp[i] - cpHOE(cpdat.temp[i], x, cpdat.natoms, cpdat.rotn));
  }
  return;
}

int cpExtrap(double *params, double *temp, double *cp, double *weight, int ndata, int natoms,
             int rotn, double tol, int maxfev, int *iwa, double *wa) {
  // iwa: len >= npar
  // wa: len >= 5*npar + npar*ndata + 2*ndata
  double ftol=tol, xtol=tol, gtol=0.0, factor=100.0, epsfcn=ndata*DBL_EPSILON;
  double *diag, *fjac, *qtf, *wa1, *wa2, *wa3, *wa4, *fvec;
  int m=ndata, n;
  int info, mode=1, nprint=0, nfev, ldfjac=m;
  int *ipvt;
  
  if (natoms == 1) { n = 3; }
  else if (natoms == 2) { n = 1; }
  else { n = 5; }
  
  ipvt = iwa;
  diag = wa; fjac = wa + n; qtf = wa + (1+m)*n;
  wa1 = wa + (2+m)*n; wa2 = wa + (3+m)*n; wa3 = wa + (4+m)*n; wa4 = wa + (5+m)*n;
  fvec = wa + (5+m)*n + m;

  cpdat.temp = temp; cpdat.cp = cp; cpdat.weight = weight;
  cpdat.ndata = ndata; cpdat.natoms = natoms; cpdat.rotn = rotn;
  
  MINPACK_LMDIF(&cpFcn, &m, &n, params, fvec, &ftol, &xtol, &gtol, &maxfev, &epsfcn, diag, &mode,
                &factor, &nprint, &info, &nfev, fjac, &ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4);
  return info;
}

double nasa7Cp(double *coeffs, double T) {
  return (coeffs[0] + coeffs[1]*T + coeffs[2]*T*T + coeffs[3]*T*T*T + coeffs[4]*T*T*T*T) * THF_RU;
}
double nasa7CpD1(double *coeffs, double T) {
  return (coeffs[1] + 2*coeffs[2]*T + 3*coeffs[3]*T*T + 4*coeffs[4]*T*T*T) * THF_RU;
}
double nasa7CpD2(double *coeffs, double T) {
  return (2*coeffs[2] + 6*coeffs[3]*T + 12*coeffs[4]*T*T) * THF_RU;
}
double nasa7H(double *coeffs, double T) {
  return (coeffs[0]*T + coeffs[1]/2.*T*T + coeffs[2]/3.*T*T*T + coeffs[3]/4.*T*T*T*T
          + coeffs[4]/5.*T*T*T*T*T + coeffs[5]) * THF_RU / 1000.;
}
double nasa7S(double *coeffs, double T) {
  return (coeffs[0]*log(T) + coeffs[1]*T + coeffs[2]/2.*T*T + coeffs[3]/3.*T*T*T
          + coeffs[4]/4.*T*T*T*T + coeffs[6]) * THF_RU;
}

int cpFit1(double *temp, double *cp, double *weight, double tempbase, double cpbase, int ndata,
           double *coeffs, double *rwork, int *iwork) {
  // iwork: len = 8*npar = 32
  // rwork: len = ndata + 2*ndata*npar + npar + npar*npar + svd_lwork <= 2342 (if ndata=231)
  double *y, *svd_a, *svd_u, *svd_s, *svd_vt, *svd_work;
  double cutoff = ndata*DBL_EPSILON;
  int *svd_iwork;
  int i, j, npar=4, svd_lwork=243, info;
  char jobz;

  y = rwork;
  svd_a = rwork + ndata;
  svd_u = rwork + ndata + ndata*npar;
  svd_s = rwork + ndata + 2*ndata*npar;
  svd_vt = rwork + ndata + 2*ndata*npar + npar;
  svd_work = rwork + ndata + 2*ndata*npar + npar + npar*npar;
  svd_iwork = iwork;
  
  for (i=0; i<ndata; i++) {
    y[i] = (cp[i] - cpbase) / THF_RU;
    svd_a[i] = weight[i] * (temp[i] - tempbase);
    svd_a[ndata+i] = weight[i] * (temp[i]*temp[i] - tempbase*tempbase);
    svd_a[2*ndata+i] = weight[i] * (temp[i]*temp[i]*temp[i] - tempbase*tempbase*tempbase);
    svd_a[3*ndata+i] = weight[i] * (temp[i]*temp[i]*temp[i]*temp[i] - tempbase*tempbase*tempbase*tempbase);
  }
  
  jobz = 'S';
  LAPACK_DGESDD(&jobz, &ndata, &npar, (double**)svd_a, &ndata, svd_s, (double**)svd_u, &ndata,
                (double**)svd_vt, &npar, svd_work, &svd_lwork, svd_iwork, &info);
  if (info != 0) { return info; }

  for (i=0; i<npar; i++) {
    svd_work[i] = 0.;
    for (j=0; j<ndata; j++) { svd_work[i] += svd_u[i*ndata+j] * y[j] * weight[j]; }
    if (svd_s[i] > cutoff) { svd_work[i] /= svd_s[i]; }
    else { svd_work[i] = 0.0; }
  }
  for (i=0; i<npar; i++) {
    coeffs[i+1] = 0.;
    for (j=0; j<npar; j++) { coeffs[i+1] += svd_work[j] * svd_vt[i*npar+j]; }
  }
  coeffs[0] = cpbase / THF_RU - (coeffs[1]*tempbase + coeffs[2]*tempbase*tempbase
                                 + coeffs[3]*tempbase*tempbase*tempbase
                                 + coeffs[4]*tempbase*tempbase*tempbase*tempbase);
  return 0;
}

int cpFit(double *temp, double *cp, int ndata, int itmin, int itmax, int itoverlap, double woverlap,
          double h298, double s298, double *coeffs, double *rwork, int *iwork) {
  // rwork: (rwork for cpFit1) + 14 + ndata <= 2587 (if ndata=231)
  double *coeffs_tmp, *weight, *rwork2;
  double mincperr=-1.0, tmid, t298, tbase, cpbase, cperr, cptmp, hmid, smid;
  double wlow=woverlap, whigh=woverlap;
  int i, j, info, it, itmid=0;
  coeffs_tmp = rwork;
  weight = rwork + 14;
  rwork2 = rwork + 14 + ndata;
  
  // find optimal tmid
  for (it=itmin; it<=itmax; it++) {
    tbase = temp[it]; cpbase = cp[it];

    // overlaped points are weighted in order to get better continuity of dCp/dT etc.
    // points near low ang high T boundaries are also weighted
    for (i=0; i<ndata; i++) {
      if ((i >= it-itoverlap) && (i <= it+itoverlap)) { weight[i] = woverlap; }
      else if (i < 2*itoverlap) { weight[i] = wlow; }
      else if (i >= ndata-2*itoverlap) { weight[i] = whigh; }
      else { weight[i] = 1.0; }
    }
    
    // fit low T points
    info = cpFit1(temp, cp, weight, tbase, cpbase, it+itoverlap, &coeffs_tmp[7], rwork2, iwork);
    if (info != 0) { return -1; }
    
    // fit high T points
    info = cpFit1(&temp[it-itoverlap], &cp[it-itoverlap], &weight[it-itoverlap], tbase, cpbase,
                  ndata-(it-itoverlap), coeffs_tmp, rwork2, iwork);
    if (info != 0) { return -1; }

    // calc sq-err
    cperr = 0.0;
    for (i=0; i<ndata; i++) {
      if (i <= it) { cptmp = nasa7Cp(&coeffs_tmp[7], temp[i]); }
      else { cptmp = nasa7Cp(coeffs_tmp, temp[i]); }
      cperr += (cptmp - cp[it]) * (cptmp - cp[it]);
    }
    if ((mincperr < 0) || (cperr < mincperr)) {
      mincperr = cperr;
      itmid = it;
      for (j=0; j<14; j++) { coeffs[j] = coeffs_tmp[j]; }
    }
  }
  
  // H & S
  t298 = 298.15;
  coeffs[12] = h298*1000./THF_RU - (coeffs[7]*t298 + coeffs[8]/2.*t298*t298
                                    + coeffs[9]/3.*t298*t298*t298
                                    + coeffs[10]/4.*t298*t298*t298*t298
                                    + coeffs[11]/5.*t298*t298*t298*t298*t298);
  coeffs[13] = s298/THF_RU - (coeffs[7]*log(t298) + coeffs[8]*t298
                              + coeffs[9]/2.*t298*t298
                              + coeffs[10]/3.*t298*t298*t298
                              + coeffs[11]/4.*t298*t298*t298*t298);
  tmid = temp[itmid];
  hmid = nasa7H(&coeffs[7], tmid);
  smid = nasa7S(&coeffs[7], tmid);
  coeffs[5] = hmid*1000./THF_RU - (coeffs[0]*tmid + coeffs[1]/2.*tmid*tmid
                                   + coeffs[2]/3.*tmid*tmid*tmid
                                   + coeffs[3]/4.*tmid*tmid*tmid*tmid
                                   + coeffs[4]/5.*tmid*tmid*tmid*tmid*tmid);
  coeffs[6] = smid/THF_RU - (coeffs[0]*log(tmid) + coeffs[1]*tmid
                             + coeffs[2]/2.*tmid*tmid
                             + coeffs[3]/3.*tmid*tmid*tmid
                             + coeffs[4]/4.*tmid*tmid*tmid*tmid);
  
  return itmid;
}

