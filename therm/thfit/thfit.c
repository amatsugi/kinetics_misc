#include "thflib.h"

int main(int argc, char *argv[]) {
  double cplst[THF_NLST], templst[THF_NLST], weightlst[THF_NLST], params[THF_NHOP];
  double cpext[THF_NEXT], tempext[THF_NEXT], coeffs[THF_NCOF];
  double *rwork;
  double econv, h298, s298, cptmp, lmtol=1e-6;
  double tmin, tmax, tmid, tnow, woverlap;
  double rmsall, rmsext, rmslst;
  double cpd1warn=0.003*THF_RU, rmswarn=0.3*THF_RU;
  int iwork[THF_IWRK];
  int maxfev=5000;
  int i, j, info, istr;
  int nlst, natoms, rotn, nvib, npar, next=THF_NEXT;
  int itmid, itmin, itmax, itoverlap;
  int badarg=0, verbose=0;
  
  FILE *ifp, *ofp;
  char ifn[THF_MAX1L], ofn[THF_MAX1L], str[THF_MAX1L], strval[THF_MAX1L];
  char name1[10+1], name2[8+1], date[6+1], elements[20+1], phase;
  char *tok;

  if (argc == 3) { strcpy(ifn, argv[1]); strcpy(ofn, argv[2]); }
  else if (argc == 4) {
    strcpy(str, argv[1]); strcpy(ifn, argv[2]); strcpy(ofn, argv[3]);
    if (strcmp(str, "-v") == 0) { verbose = 1; }
    else if (strcmp(str, "--verbose") == 0) { verbose = 1; }
    else { printf("Unknown option: %s\n", str); badarg = 1; }
  }
  else { badarg = 1; }

  if (badarg > 0)
    { printf("Usage: thfit [-v|--verbose] input_file output_file\n"); return 1; }

  if ((ifp = fopen(ifn, "r")) == NULL)
    { printf("Error: failed to open %s\n", ifn); return 1; }
  if ((ofp = fopen(ofn, "w")) == NULL)
    { printf("Error: failed to open %s\n", ofn); fclose(ifp); return 1; }
  
  rwork = (double*)malloc(sizeof(double)*THF_RWRK);
  if (rwork == NULL)
    { printf("failed to allocate array\n"); fclose(ifp); fclose(ofp); return 1; }
  
  // prep. tempext
  tmin = 300.; tmax = 5000.;
  itmin = 50; itmax = 120; // 800 <= tmid <= 1500
  itoverlap = 5; woverlap = 5.; // +- 5 overlap points (weighted by 5)
  for (i=0; i<next; i++) {
    if (i < 170) { tempext[i] = 300. + 10. * i; } // 300-2000, step=10
    else { tempext[i] = 2000. + 50. * (i-170); } // 2000-5000, step=50
  }
  
  // templst
  templst[0] =  300.0; weightlst[0] = 1.0;
  templst[1] =  400.0; weightlst[1] = 1.0;
  templst[2] =  500.0; weightlst[2] = 1.0;
  templst[3] =  600.0; weightlst[3] = 1.0;
  templst[4] =  800.0; weightlst[4] = 1.0;
  templst[5] = 1000.0; weightlst[5] = 1.0;
  templst[6] = 1500.0; weightlst[6] = 1.0;

  // header
  fprintf(ofp, "THERMO\n");
  fprintf(ofp, "   300.000  1500.000  5000.000\n");
  
  // first line: unit
  if ((fgets(str, THF_MAX1L-1, ifp)) == NULL)
    { printf("Error: failed to read UNITS\n"); fclose(ifp); fclose(ofp); free(rwork); return 1; }
  
  istr = 0;
  for (i=0; i < (int)strlen(str); i++) { if (str[i] != ' ') { istr = i; break; } }
  for (i=istr; i<(int)strlen(str); i++) {
    if ((str[i] == ' ') || (str[i] == '\t') || (str[i] == '\r') || (str[i] == '\n'))
      { str[i] = '\0'; break; }
  }
  if (strcmp(&str[istr], "UNITS:KJ") == 0) { econv = 1.0; }
  else if (strcmp(&str[istr], "UNITS:KCAL") == 0) { econv = THF_CAL2J; }
  else { printf("Error: Unknown UNITS: %s\n", &str[istr]); fclose(ifp); fclose(ofp); free(rwork); return 1; }
  
  // second line: title
  if ((fgets(str, THF_MAX1L-1, ifp)) == NULL)
    { printf("Error: reading title\n"); fclose(ifp); fclose(ofp); free(rwork); return 1; }
  // third line: header
  if ((fgets(str, THF_MAX1L-1, ifp)) == NULL)
    { printf("Error: reading header\n"); fclose(ifp); fclose(ofp); free(rwork); return 1; }
  
  // lst
  while ((fgets(str, THF_MAX1L-1, ifp)) != NULL) {
    if ((int)strlen(str) < 131)
      { printf("Error: reading data\n"); fclose(ifp); fclose(ofp); free(rwork); return 1; }
    // read line
    str[130] = '\0';
    strncpy(name1, &str[1], 10); name1[10] = '\0';
    strncpy(name2, &str[86], 8); name2[8] = '\0';
    strncpy(date, &str[94], 6); date[6] = '\0';
    strncpy(elements, &str[103], 5); strncpy(&elements[5], &str[109], 5);
    strncpy(&elements[10], &str[115], 5); strncpy(&elements[15], &str[121], 5);
    elements[20] = '\0';
    phase = str[127];
    rotn = atoi(&str[128]);

    // h, s, cp
    strncpy(strval, &str[12], 73); strval[73] = '\0';
    h298 = atof(strtok(strval, " ")) * econv;
    s298 = atof(strtok(NULL, " ")) * econv;
    for (i=0; i<6; i++) { cplst[i] = atof(strtok(NULL, " ")) * econv; }
    tok = strtok(NULL, " ");
    if (tok == NULL) { nlst = 6; }
    else { cplst[6] = atof(tok) * econv; nlst = 7; }

    // calc natoms
    natoms = 0;
    strncpy(strval, &elements[2], 3); strval[3] = '\0'; tok = strtok(strval, " ");
    if (tok != NULL) { natoms += atoi(tok); }
    strncpy(strval, &elements[7], 3); strval[3] = '\0'; tok = strtok(strval, " ");
    if (tok != NULL) { natoms += atoi(tok); }
    strncpy(strval, &elements[12], 3); strval[3] = '\0'; tok = strtok(strval, " ");
    if (tok != NULL) { natoms += atoi(tok); }
    strncpy(strval, &elements[17], 3); strval[3] = '\0'; tok = strtok(strval, " ");
    if (tok != NULL) { natoms += atoi(tok); }
  
    printf("Processing %s%s\n", name1, name2);
    if (verbose > 0) {
      printf("  date: %s\n", date);
      printf("  rotn: %d\n", rotn);
      printf("  phase: %c\n", phase);
      printf("  elements: %s\n", elements);
      printf("  natoms: %d\n", natoms);
    }
    
    // guess
    if (natoms == 1) { npar = 3; params[0] = 0.8; params[1] = 500.; params[2] = 1000.; }
    else if (natoms == 2) { npar = 1; params[0] = 1000.; }
    else {
      npar = 5;
      nvib = 3 * natoms - 6 - rotn;
      params[0] = nvib*0.3; params[1] = nvib*0.4;
      params[2] = 800.; params[3] = 1600.; params[4] = 3000.;
    }
    
    // extrapolate
    if ((natoms == 1) && (fabs(cplst[0] - cplst[nlst-1]) < 0.01)) {
      // constant cp
      for (i=0; i<next; i++) { cpext[i] = cplst[0]; }
      printf("Extrapolate: constant Cp\n");
    }
    else {
      // fit HOE
      printf("Extrapolating by HOE...\n");
      info = cpExtrap(params, templst, cplst, weightlst, nlst, natoms,
                      rotn, lmtol, maxfev, iwork, rwork);
      if ((info == 0) || (info >= 5))
        { printf("Error: LMDIF info = %d\n", info); fclose(ifp); fclose(ofp); free(rwork); return 1; }
      if (verbose > 0) {
        printf("  LMDIF info = %d\n", info);
        printf("    Fitted params:");
        for (i=0; i<npar; i++) { printf(" %g", params[i]); }
        printf("\n");
      }
      for (i=0; i<next; i++) { cpext[i] = cpHOE(tempext[i], params, natoms, rotn); }
    }
    
    // fit NASA7 format
    printf("Fitting to NASA7 format...\n");
    itmid = cpFit(tempext, cpext, next, itmin, itmax, itoverlap, woverlap,
                  h298, s298, coeffs, rwork, iwork);
    if (itmid < 0)
      { printf("Error: DGESDD did not converged.\n"); fclose(ifp); fclose(ofp); free(rwork); return 1; }
    tmid = tempext[itmid];

    // check continuity of dCp/dT
    if (fabs(nasa7CpD1(&coeffs[7], tmid) - nasa7CpD1(coeffs, tmid)) > cpd1warn) {
      printf("Warning: large discontinuity of dCp/dT at Tmid:\n");
      printf("  Tmid = %10.3f\n", tmid);
      printf("  dCp/dT = %12.6f %12.6f\n", nasa7CpD1(&coeffs[7], tmid), nasa7CpD1(coeffs, tmid));
    }
    if (verbose > 0) {
      printf("  Tmid = %10.3f\n", tmid);
      printf("  At Tmid:       low          high\n");
      printf("    Cp       = %12.6f %12.6f\n", nasa7Cp(&coeffs[7], tmid), nasa7Cp(coeffs, tmid));
      printf("    dCp/dT   = %12.6f %12.6f\n", nasa7CpD1(&coeffs[7], tmid), nasa7CpD1(coeffs, tmid));
      printf("    d2Cp/dT2 = %12.6f %12.6f\n", nasa7CpD2(&coeffs[7], tmid), nasa7CpD2(coeffs, tmid));
    }

    // check extrap. & fitting error
    rmsall = 0.; rmsext = 0.; rmslst = 0.;
    if (verbose > 0) { printf("     T[K]      Cp(ext)      Cp(fit)      Cp(lst) [JK-1mol-1]\n"); }
    tnow = 300.;
    for (i=0; i<next; i++) {
      if (i <= itmid) { cptmp = nasa7Cp(&coeffs[7], tempext[i]); }
      else { cptmp = nasa7Cp(coeffs, tempext[i]); }
      rmsall += (cptmp - cpext[i])*(cptmp - cpext[i]);
      if (fabs(tnow - tempext[i]) < 0.001) {
        if (verbose > 0) { printf(" %10.3f %12.6f %12.6f", tempext[i], cpext[i], cptmp); }
        for (j=0; j<nlst; j++) {
          if (fabs(templst[j] - tempext[i]) < 0.001) {
            rmsext += (cpext[i] - cplst[j])*(cpext[i] - cplst[j]);
            rmslst += (cptmp - cplst[j])*(cptmp - cplst[j]);
            if (verbose > 0) { printf(" %12.6f", cplst[j]); }
          }
        }
        if (verbose > 0) { printf("\n"); }
        if (tnow > 1999.) { tnow += 500; }
        else { tnow += 100.; }
      }
    }
    rmsall = sqrt(rmsall / (float)next);
    rmsext = sqrt(rmsext / (float)nlst);
    rmslst = sqrt(rmslst / (float)nlst);
    printf("RMS fitting errors:\n");
    printf("  Extrap.  (n=%3d): %10.6f [JK-1mol-1]\n", nlst, rmsext);
    printf("  FitNASA7 (n=%3d): %10.6f [JK-1mol-1]\n", next, rmsall);
    printf("  Ext&Fit  (n=%3d): %10.6f [JK-1mol-1]\n", nlst, rmslst);
    if ((rmsall > rmswarn) || (rmsext > rmswarn) || (rmslst > rmswarn))
      { printf("Warning: large RMS error\n"); }
    printf("\n");
    
    // output
    fprintf(ofp, "%s%s%s%s%c", name1, name2, date, elements, phase);
    fprintf(ofp, "% 10.3f% 10.3f% 8.2f      1\n", tmin, tmax, tmid);
    for (i=0; i<=4; i++) { fprintf(ofp, "% 15.7E", coeffs[i]); }
    fprintf(ofp, "    2\n");
    for (i=5; i<=9; i++) { fprintf(ofp, "% 15.7E", coeffs[i]); }
    fprintf(ofp, "    3\n");
    for (i=10; i<=13; i++) { fprintf(ofp, "% 15.7E", coeffs[i]); }
    fprintf(ofp, "                   4\n");
  }
  
  // footer
  fprintf(ofp, "END\n");
  
  fclose(ofp);
  fclose(ifp);
  free(rwork);
  
  return 0;
}

