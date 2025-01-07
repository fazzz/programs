
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "WHAM.h"
#include "REMD.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

int main(int argc, char *argv[]) {
  int i,j,k;

  int nwindows,numatom,*nt;
  int nttemp;
  double *ebF,***ebW;
  int **dihedpairs;
  int numdrest;
  double **theta_ref,V_rest;
  double ***trj;
  double  temp;

  char	 *inputfilename1,*inputfilename2,*parmfilename;
  char	 *outputfilename,*outputfilename2;

  FILE	*inputfile1,*inputfile2,*outputfile,*outputfile2,*parmfile;

  nwindows=atoi(getenv("nwindows"));
  nt=(int *)gcemalloc(sizeof(int)*nwindows);
  for (i=0;i<nwindows;++i) {
    sprintf(nti,"nt%d",i+1);
    nt[i]=atoi(getenv(nti));
  }
  for (i=0;i<nwindows;++i) {
    sprintf(topi,"top%d",i+1);
    parmfiles[i]=atoi(getenv(topi));
  }
  parmfilename=getenv("top_ref");
  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;
  outputfilename=*++argv;
  trj=(double ***)gcemalloc(sizeof(double **)*nwindows);
  for (i=0;i<nwindows;++i) {
    trj[i]=(double **)gcemalloc(sizeof(double *)*nt[i]);
    for (j=0;j<nt[i];++j)
      trj[i][j]=(double *)gcemalloc(sizeof(double)*numatom*3);
  }
  for (i=0;i<nwindows;++i) {
    sprintf(trji,"trj%d",i+1);
    inputfilename1=atoi(getenv(trji));
    inputfile1=efopen(inputfilename1,"r");
    io_scantrj(inputfile1,numatom,nt[i],trj[i]);
    fclose(inputfile1);
  }

  ebF  = (double *)gcemalloc(sizeof(double)*nwindows);
    
  ebW  = (double ***)gcemalloc(sizeof(double **)*nwindows);
  for (i=0;i<nwindows;++i) {
    ebW[i]=(double **)gcemalloc(sizeof(double *)*nwindows);
    for (j=0;j<nwindows;++j) {
      ebW[i][j]=(double *)gcemalloc(sizeof(double)*nt[j]);
    }
  }
  
  wham_calc_force_SDFF(nwindows,nt,ebW,trj,numatom,parmtops,parmtop_ref,temp);
  wham_ite(nwindows,nt,ebF,ebW);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<nwindows;++i)
    fprintf(outputfile,"%e\n",ebF[i]);  
  fclose(outputfile);

  outputfile2=efopen("log_ef.txt","w");
  for (i=0;i<nwindows;++i) {
    for (j=0;j<nwindows;++j) {
      for (k=0;k<nt[j];++k)
	fprintf(outputfile2,"%e ",ebW[i][j][k]);
      fprintf(outputfile2,"\n");
    }
  }
  fclose(outputfile2);
 
  return 0;
}

