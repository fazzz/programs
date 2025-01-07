
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EMalg.h"
#include "K_means.h"
#include "Gaussian.h"

#include "STRING.h"
#include "CSI.h"
#include "PT.h"
#include "FF.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

int USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;

  int sl=0;

  int outinterval=1;

  int numpoint,numiteration=10000;
  int numCV;

  double dt=0.01;

  double maxdelta,*pe_old;
  double criteria=1.0e-3;

  char *line;
  size_t len=0;

  char *progname;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *path,*path_evoluted;
  double *fe,*pe;

  numpoint=atoi(*++argv);
  inifilename=*++argv;

  path=(double *)gcemalloc(sizeof(double)*numCV*numpoint);
  inifile=efopen(inifilename,"r");
  for (i=0;i<sl;++i) getline(&line,&len,inifile);
  io_scantrj2(inifile,numatom,numpoint,path);
  fclose(inifile);

  path_evoluted=(double *)gcemalloc(sizeof(double)*numCV*numpoint);
  fe=(double *)gcemalloc(sizeof(double)*numCV*numpoint);
  pe=(double *)gcemalloc(sizeof(double)*numpoint);
  pe_old=(double *)gcemalloc(sizeof(double)*numpoint);
  crd=(double *)gcemalloc(sizeof(double)*numCV);

  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");
  logfile1=efopen("log1_zswoH","w");
  logfile2=efopen("log2_zswoH","w");
  for (i=0;i<numiteration;++i) {
    for (j=0;j<numpoint*numatom*3;++j) fe[j]=0.0;
    for (j=0;j<numpoint;++j) pe[j]=0.0;

    for (j=0;j<numpoint;++j) {
      for (k=0;k<numCV;++k) {
	crd[k]=path[j*CV+k];

	de_twoD_Mixed_Gaussian(crd,f,nyu,Sigma,pi);

    }
    z_string_Cartesian(path,path_evoluted,fe,numpoint,numCV,dt);

    printf("%d %lf %lf %lf\n",i,pe[0],pe[(int)(numpoint/2)],pe[numpoint-1]);

    if (i>0) {
      maxdelta=fabs(pe[0]-pe_old[0]);
      for (j=1;j<numpoint;++j) {
	if (maxdelta<abs(pe[j]-pe_old[j])) maxdelta=pe[j]-pe_old[j];
      }
      printf("%lf\n",fabs(maxdelta));
      if (fabs(maxdelta)<criteria) break;
    }
    for (j=0;j<numpoint;++j) pe_old[j]=pe[j];
    if (i%outinterval==0) {
      io_outtrj2(logfile1,numatom,numpoint,path);
      for (j=0;j<numpoint;++j) fprintf(logfile2,"%4d %12.8e \n",j,pe[j]);
    }
  }
  fclose(logfile1);
  fclose(logfile2);

  io_outtrj2(outputfile,numatom,numpoint,path);
  fprintf(outputfile2,"# ene\n");
  for (i=0;i<numpoint;++i) fprintf(outputfile2,"%4d %12.8lf\n",i,pe[i]);
  fclose(outputfile);
  fclose(outputfile2);

}

int USAGE(char *progname) {
  printf("USAGE:%s \n",progname);
  printf("[-t dt]    \n");
  printf("[-o outinterval] \n");
  printf("[-c criteria] \n");
  printf("[-m numiteration] \n");
  printf("[-h] help  \n");
  printf("%s numpoint inifilename(path) parmtopname outputfilename(path) outputfilename2(ene)\n",progname);
}
