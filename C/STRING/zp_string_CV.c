
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "STRING_CV.h"
#include "STRING.h"
#include "EF.h"

void z_string_CV(double **path,double **path_evoluted,double **fe,int numpoint,int numCV,double dt) {
  int i,j,k;
  double *path2,*path_evoluted2,*fe2;

  path2=(double *)gcemalloc(sizeof(double)*numpoint*numCV);
  path_evoluted2=(double *)gcemalloc(sizeof(double)*numpoint*numCV);
  fe2=(double *)gcemalloc(sizeof(double)*numpoint*numCV);

  k=0;
  for (i=0;i<numpoint;++i) {
    for (j=0;j<numCV;++j) {
      path2[k]=path[i][j];
      path_evoluted2[k]=path_evoluted[i][j];
      fe2[k]=fe[i][j];
      ++k;
    }
  }

  z_string_evolution(path2,path_evoluted2,fe2,numpoint,dt,numCV);
  z_string_interpolation(path2,path_evoluted2,numpoint,numCV);

  k=0;
  for (i=0;i<numpoint;++i) {
    for (j=0;j<numCV;++j) {
      path[i][j]=path2[k];
      path_evoluted[i][j]=path_evoluted2[k];
      ++k;
    }
  }

}
