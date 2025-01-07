
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "PT.h"
#include "dPCA.h"
#include "IO.h"
#include "EF.h"

int main(int argc, char *argv[]) {
  int i,j,c;
  int flag=0;
  int numstep,numdihed;
  double sum,sum2;
  double *dtrj,*sctrj,*cov,*eigenval,*dpca,tp;

  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename1,*parmfilename;
  char *outputfilename,*outputfilename2,*outputfilename3;
  FILE *inputfile1,*parmfile;
  FILE *outputfile,*outputfile2,*outputfile3,*logfile;

  /*******************************************/
  /* while((c=getopt(argc,argv,"FL"))!=-1) { */
  /*   switch(c) {			     */
  /*   case 'F':			     */
  /*     flag='c';			     */
  /*     break;				     */
  /*   case 'L':			     */
  /*     flag='x';			     */
  /*     break;				     */
  /*   default:				     */
  /*     usage();			     */
  /*     exit(1);			     */
  /*   }				     */
  /* }					     */
  /* 					     */
  /* argc -= optind;			     */
  /* argv += optind;			     */
  /*******************************************/
  
  if (argc < 8) {
    printf("USAGE: %s numdihed tp numstep inputfilename1(data) inputfilename2(pt) outputfilename(log) outputfilename2(mode) outputfilename3(dpca)\n",argv[0]);
    exit(1);
  }
  numdihed=atoi(*++argv);
  tp=atof(*++argv);
  numstep=atoi(*++argv);
  inputfilename1 = *++argv;
  parmfilename = *++argv;
  outputfilename = *++argv;
  outputfilename2= *++argv;
  outputfilename3= *++argv;
  
  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);

  dtrj=(double *)gcemalloc(sizeof(double)*numstep*numdihed);
  sctrj=(double *)gcemalloc(sizeof(double)*numstep*numdihed*2);
  dpca=(double *)gcemalloc(sizeof(double)*numstep*numdihed*2);
  cov=(double *)gcemalloc(sizeof(double)*numdihed*2*numdihed*2);
  eigenval=(double *)gcemalloc(sizeof(double)*numdihed*2);

  for (i=0;i<numdihed*2;++i)
    for (j=0;j<numdihed*2;++j)
      cov[i*numdihed*2+j]=0.0;
  inputfile1=efopen(inputfilename1,"r");
  io_scandtraj(inputfile1,numstep,numdihed,dtrj);
  fclose(inputfile1);

  dpca_norm(dtrj,sctrj,numstep,numdihed);
  dpca_covm(sctrj,numstep,numdihed,cov);
  dpca_diag(cov,eigenval,numdihed);
  dpca_proj(sctrj,dpca,cov,numstep,numdihed);
  
  sum=0.0;sum2=0;
  for (i=0;i<numdihed*2;++i)
    sum+=eigenval[i];
  outputfile=efopen(outputfilename,"w");
  fprintf(outputfile,"n eigenvalue(sigma^2) %% \n");
  for (i=0;i<numdihed*2;++i) {
    sum2+=eigenval[i];
    fprintf(outputfile,"%d %12.8lf %12.8lf \n",i,eigenval[i],sum2/sum*100.0);
  }
  fclose(outputfile);

  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<numdihed*2;++i) {
    for (j=0;j<numdihed*2;++j)
      fprintf(outputfile2,"%lf ",cov[i*numdihed*2+j]);
    fprintf(outputfile2,"\n");
  }
  fclose(outputfile2);

  outputfile3=efopen(outputfilename3,"w");
  for (i=0;i<numstep;++i) {
    fprintf(outputfile3,"%d ",i);
    for (j=0;j<numdihed*2;++j)
      fprintf(outputfile3,"%lf ",dpca[i*numdihed*2+j]);
    fprintf(outputfile3,"\n");
  }
  fclose(outputfile3);

  logfile=efopen("log_sc.txt","w");
  for (i=0;i<numstep;++i) {
    fprintf(logfile,"%d ",i);
    for (j=0;j<numdihed*2;++j)
      fprintf(logfile,"%lf ",sctrj[i*numdihed*2+j]);
    fprintf(logfile,"\n");
  }
  fclose(logfile);

  return 0;
}
