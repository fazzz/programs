
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <math.h>

#include "EF.h"

#include "Yacobian.h"

#define ON 0
#define OFF 1

double Yacobian_wrtdihed(double *crd,struct AmberParmL ap,
			 int numatom,int na1,int na2,double *Yac) {
  int i,j,k;
  int MSflag;
  int LastAtom;

  double l,*e;

  e=(double *)gcemalloc(sizeof(double)*3);

  for (i=0;i<3;++i) e[i]=crd[na2*3+i]-crd[na1*3+i];

  l=0.0; for (i=0;i<3;++i) l+=e[i]*e[i]; l=sqrt(l);

  for (i=0;i<3;++i) e[i]=e[i]/l;

  LastAtom=numatom;

  if (strcmp(ap.ITREE[na2],"M")==0) {
    MSflag=ON;
  }
  else  { 
    MSflag=OFF;
    for (i=na2;i<numatom;++i) {
      if (strcmp(ap.ITREE[i],"E")==0) {
	LastAtom=i;
	break;
      }
    }
  }

  for (i=0;i<numatom;++i) {
    if (i>na2 && i<= LastAtom) {
      Yac[i*3]=crd[i*3+2]*e[1]-crd[i*3+1]*e[2];
      Yac[i*3+1]=crd[i*3]*e[2]-crd[i*3+2]*e[0];
      Yac[i*3+2]=crd[i*3+1]*e[0]-crd[i*3]*e[1];
    }
    else {
      for (j=0;j<3;++j) Yac[i*3+j]=0.0;
    }
  }
  
  return 1.0;
}
