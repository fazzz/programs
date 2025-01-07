
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "PT.h"
#include "EF.h"

#define UNIT 4.184070*100.0

double kieneftrj(double *crd1, double *crd2, double *mass, double dt, int numatom) {
  int i,j;
  double ke=0.0;

  for (i=0;i<numatom*3;++i) {
    j=(int)(i/3);
    ke+=0.5*mass[j]*(crd1[i]-crd2[i])*(crd1[i]-crd2[i])/(dt*dt);
  }

  return ke/UNIT;
}

