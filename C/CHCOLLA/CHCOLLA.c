
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PT.h"
#include "CB.h"
#include "CA.h"
#include "mymath.h"
#include "TOPO.h"

double CHCO_diff_bond_lengh(double *crd, double *crdref, int MODE, double *ave, double *var) {
  int i;
  double *bl,*blref,*dbl;
  int numb;

  if (MODE==INCH){  
    bl=(double *)gcemalloc(sizeof(double)*(AP.MBONA+AP.NBONH));
    blref=(double *)gcemalloc(sizeof(double)*(AP.MBONA+AP.NBONH));
    dbl=(double *)gcemalloc(sizeof(double)*(AP.MBONA+AP.NBONH));
  }
  else if (MODE==EXCH) {
    bl=(double *)gcemalloc(sizeof(double)*AP.MBONA);
    blref=(double *)gcemalloc(sizeof(double)*AP.MBONA);
    dbl=(double *)gcemalloc(sizeof(double)*AP.MBONA);
  }
  numb=CB(crd,MODE,bl);
  CB(crdref,MODE,blref);

  for (i=0;i<numb;++i) dbl[i]=fabs(bl[i]-blref[i]);

  *ave=calc_ave(numb,dbl);
  *var=calc_var(numb,dbl);

}

double CHCO_diff_bond_angle(double *crd, double *crdref, int MODE, double *ave, double *var) {
  int i;
  double *ba,*baref,*dba;
  double pi;
  int numa;

  pi=acos(-1.0);
  if (MODE==INCH){  
    ba=(double *)gcemalloc(sizeof(double)*(AP.MTHETA+AP.NTHETH));
    baref=(double *)gcemalloc(sizeof(double)*(AP.MTHETA+AP.NTHETH));
    dba=(double *)gcemalloc(sizeof(double)*(AP.MTHETA+AP.NTHETH));
  }
  else if (MODE==EXCH) {
    ba=(double *)gcemalloc(sizeof(double)*AP.MTHETA);
    baref=(double *)gcemalloc(sizeof(double)*AP.MTHETA);
    dba=(double *)gcemalloc(sizeof(double)*AP.MTHETA);
  }
  numa=CA(crd,MODE,ba);
  CA(crdref,MODE,baref);

  for (i=0;i<numa;++i) {
    dba[i]=ba[i]-baref[i];
    if (dba[i] < 2.0*pi-dba[i])
      dba[i]=2.0*pi-dba[i];
  }

  *ave=calc_ave(numa,dba);
  *var=calc_var(numa,dba);

}

