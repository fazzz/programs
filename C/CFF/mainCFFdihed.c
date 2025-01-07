
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "FF.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

int main(int argc, char *argv[]) {
  int i,j,k,flag;
  int f;

  double *p_d;
  double *n_d;
  int numatom,numstep;
  double *crd;

  char *inputfilename,*inputfilename2;
  char *outputfilename;
  FILE *inputfile,*inputfile2;
  FILE *outputfile;
  FILE *log;

  if (argc < 4) {
    printf("USAGE: %s numstep inputfilename(crd) inputfilename2(top) outputfilename(ene_di) \n",argv[0]);
    exit(1);
  }
  numstep=atoi(*++argv);
  inputfilename   = *++argv;
  inputfilename2  = *++argv;
  outputfilename  = *++argv;
  
  inputfile2=efopen(inputfilename2,"r");
  readParmtop(inputfile2);
  fclose(inputfile2);
  numatom=AP.NATOM;

  p_d=(double *)gcemalloc(sizeof(double)*(AP.NPHIH+AP.MPHIA));
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfile=efopen(inputfilename,"r");
  outputfile=efopen(outputfilename,"w");

  for (i=0;i<numstep;++i) {
    io_scanconf(inputfile,numatom,crd,'x');    
    ff_calcDIHE(p_d,n_d,crd,1,0,0);

    fprintf(outputfile,"%d ",i);
    for (j=0;j<AP.NPHIH;++j)
      fprintf(outputfile,"%e ",p_d[j]);
    for (j=0;j<AP.MPHIA;++j)
      fprintf(outputfile,"%e ",p_d[j+AP.NPHIH]);
      fprintf(outputfile,"\n");
  }

  log=efopen("log_CFFDI.txt","w");
  for (i=0;i<AP.NPHIH;++i) {
    fprintf(log,"%d %d %d %d %d\n",abs(AP.PH[i][0])/3+1,abs(AP.PH[i][1])/3+1,abs(AP.PH[i][2])/3+1,abs(AP.PH[i][3])/3+1,abs(AP.PH[i][4]));
  }
  for (i=0;i<AP.MPHIA;++i) {
    fprintf(log,"%d %d %d %d %d\n",abs(AP.PA[i][0])/3+1,abs(AP.PA[i][1])/3+1,abs(AP.PA[i][2])/3+1,abs(AP.PA[i][3])/3+1,abs(AP.PA[i][4]));
  }

  fclose(log);
  fclose(inputfile);
  fclose(outputfile);

  return 0;
}

