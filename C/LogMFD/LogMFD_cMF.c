
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <math.h>

#include "EF.h"

#include "Yacobian.h"

#include "LogMFD_cMF.h"


double LogMFD_cdF_dX_wrtdihed(double *crd,struct potential *e, struct force *f, struct AmberParmL ap,
			      int numatom, int na1, int na2) {
  int i,j,k;

  double *ff,*Yac;
  double dF_dX;

  ff=(double *)gcemalloc(sizeof(double)*numatom*3);
  Yac=(double *)gcemalloc(sizeof(double)*numatom*3);

  ffLc_calcffandforce(crd,numatom,e,f,ap);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      ff[i*3+j]=+(*f).f_b[i*3+j]-(*f).f_a[i*3+j]-(*f).f_d[i*3+j]
	-(*f).f_e[i*3+j]-(*f).f_LJ[i*3+j]-(*f).f_e_14[i*3+j]-(*f).f_LJ_14[i*3+j];

  Yacobian_wrtdihed(crd,ap,numatom,na1,na2,Yac);

  dF_dX=0.0; for (i=0;i<numatom*3;++i) dF_dX+=Yac[i]*ff[i];

  return dF_dX;
}


