
#include <stdio.h>

#include "RADG.h"

double RADG_calc_radg(double **crd, 
		      double *mass,
		      int numatom ) {
  int i,j;
  double com[3],smass=0.0,RADG=0.0;

  for (i=0;i<3;++i) {
    com[i]=0.0;
  }

  for (i=0;i<numatom;++i) {
    smass+=mass[i];
  }
  
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      com[j]+=mass[i]*crd[i][j];
    }
  }

  for (i=0;i<3;++i) {
    com[i]=com[i]/smass;
  }

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      RADG+=sqrt(mass[i]*(crd[i][j]-com[j])*(crd[i][j]-com[j])/smass);
    }
  }

  return RADG;
}
