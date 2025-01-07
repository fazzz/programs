
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "WHAM.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

int main(int argc, char *argv[]) {
  int	i,j,k;

  int		  nwindows,numatom,*nt,*dummy,numnb;
  double	 *ebF,***ebW,***trj;
  double	**ui;
  double	 *fact;
  double	  temp;

  char	  inputfilename[100],inputfilename2[100],*parmfilename;
  char	 *outputfilename,*outputfilename2;
  char	**uifilename,**trjinputfilename,nti[100];

  FILE	*inputfile,*inputfile2,*outputfile,*outputfile2,*parmfile;

  nwindows=atoi(getenv("nwindows"));
  nt   = (int *)gcemalloc(sizeof(int)*nwindows);
  dummy= (int *)gcemalloc(sizeof(int)*nwindows);
  parmfilename=getenv("parmfilename");
  outputfilename=getenv("outputfilename");
  outputfilename2=getenv("outputfilename2");
  uifilename = (char **)gcemalloc(sizeof(char *)*nwindows);
  trjinputfilename = (char **)gcemalloc(sizeof(char *)*nwindows);
  for (i=0;i<nwindows;++i) {
    sprintf(inputfilename,"inputfilename1%d",i+1);
    sprintf(inputfilename2,"inputfilename2%d",i+1);
    sprintf(nti,"nt%d",i+1);
    trjinputfilename[i]=getenv(inputfilename);
    uifilename[i]=getenv(inputfilename2);
    dummy[i]=atoi(getenv(nti));
  }
  memcpy(nt,dummy,sizeof(int)*nwindows);
  temp=atoi(getenv("temp"));

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;
  numnb=ff_set_numnb();

  fact = (double *)gcemalloc(sizeof(int)*nwindows);
  ui=(double **)gcemalloc(sizeof(double *)*nwindows);
  for (i=0;i<nwindows;++i)
    ui[i]=(double *)gcemalloc(sizeof(double)*numnb*2);
  for (i=0;i<nwindows;++i) {
    inputfile=efopen(uifilename[i],"r");
    for (j=0;j<numnb;++j)
      fscanf(inputfile,"%lf",&ui[i][j*2]);
    for (j=0;j<numnb;++j)
      fscanf(inputfile,"%lf",&ui[i][j*2+1]);
    fscanf(inputfile,"%lf",&fact[i]);
    fclose(inputfile);
  }

  ebF  = (double *)gcemalloc(sizeof(double)*nwindows);

  trj  = (double ***)gcemalloc(sizeof(double **)*nwindows);
  for (i=0;i<nwindows;++i) {
    trj[i]=(double **)gcemalloc(sizeof(double *)*nt[i]);
    for (j=0;j<nt[i];++j) {
      trj[i][j]=(double *)gcemalloc(sizeof(double)*numatom*3);
    }
  }
    
  ebW  = (double ***)gcemalloc(sizeof(double **)*nwindows);
  for (i=0;i<nwindows;++i) {
    ebW[i]=(double **)gcemalloc(sizeof(double *)*nwindows);
    for (j=0;j<nwindows;++j) {
      ebW[i][j]=(double *)gcemalloc(sizeof(double)*nt[j]);
    }
  }

  for (i=0;i<nwindows;++i) {
    inputfile=efopen(trjinputfilename[i],"r");
    io_scantrj(inputfile,numatom,nt[i],trj[i]);
    fclose(inputfile);
  }
  
  wham_calc_force(nwindows,nt,ebW,trj,ui,fact,temp);
  wham_ite(nwindows,nt,ebF,ebW);
  //  wham_pmf(num_k,num_x,num_y,num_window,nt,ebF,ebW,hist,pmf);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<nwindows;++i)
    fprintf(outputfile,"%d %e\n",i,ebF[i]);  
  fclose(outputfile);
 
  return 0;
}

