
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "EF.h"
#include "PT.h"
#include "CFF.h"

void simplex_reflection(double *parametersets/**/,double *newparameterset, int numofparamaters,int indexofmax) {
  int i,j;

  for (i=0;i<numofparamaters;++i) {
    newparameterset[i] =0.0;
  }

  for (i=0;i<numofparamaters;++i) {
    for (j=0;j<numofparamaters+1;++j) {
      newparameterset[i] += 2.0/numofparamaters*parametersets[j];
    }
    newparameterset[i]-=(2.0/numofparamaters+1.0)*parametersets[indexofmax];
  }

}




