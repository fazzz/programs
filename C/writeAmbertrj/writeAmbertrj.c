#include <stdio.h>
#include <stdlib.h>

#include "EF.h"

#include "writeAmbertrj.h"

int writeAmbertrj(double *crd, int numatom, FILE *Ambertrj) {  
  int i,j,k,l;

  k=1;
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      fprintf(Ambertrj," %6.3lf",crd[i*3+j]);
      if ((l=(k%10))==0) fprintf(Ambertrj,"\n");
      ++k;
    }
  }

}

int makeAmbertrj(char *Ambertrjfilename,char *protname, FILE *Ambertrj) {  

  Ambertrj=efopen(Ambertrjfilename,"a");
  
  fprintf(Ambertrj," %4s\n",protname);

}
