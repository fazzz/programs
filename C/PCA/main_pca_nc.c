
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PCA.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"
#include "netcdf_mine.h"
//#include "netcdf_mineL.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,s;
  int spdim=2;
  int numatom,numatomp,numstep;
  double *traj,*cov,*eigenvalue;

  int MODE=AA,/*inMODE=AA,outMODE=CA,*/IOMODE=MD;

  char *inputfilename,*parmtopfilename;
  char *outputfilename1,*outputfilename2,*outputfilename3;
  char *progname;

  FILE *inputfile, *parmtop;
  FILE *outputfile1,*outputfile2,*outputfile3;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  int opt_idx=1;

  struct option long_opt[] = {
    /*************************/
    /* {"inCA",0,NULL,'c'},  */
    /* {"inAA",0,NULL,'a'},  */
    /* {"ineH",0,NULL,'e'},  */
    /* {"outCA",0,NULL,'C'}, */
    /* {"outAA",0,NULL,'A'}, */
    /* {"outeH",0,NULL,'E'}, */
    /*************************/
    {"b",0,NULL,'b'},
    {"c",0,NULL,'c'},
    {"a",0,NULL,'a'},
    {"K",0,NULL,'K'},
    {"s",1,NULL,'s'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"bcaKhs:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'K':
      IOMODE=AMBER;
      break;
    /*****************/
    /* case 'c':     */
    /*   inMODE=CA;  */
    /*   break;	     */
    /* case 'a':     */
    /*   inMODE=AA;  */
    /*   break;	     */
    /* case 'e':     */
    /*   inMODE=HV;  */
    /*   break;	     */
    /* case 'C':     */
    /*   outMODE=CA; */
    /*   break;	     */
    /* case 'A':     */
    /*   outMODE=AA; */
    /*   break;	     */
    /* case 'E':     */
    /*   outMODE=HV; */
    /*   break;	     */
    /*****************/
    case 'b':
      MODE=CA;
      break;
    case 'c':
      MODE=HV;
      break;
    case 'a':
      MODE=AA;
      break;
    case 's':
      spdim=atoi(optarg);
      break;
    case 'h':
      USAGE(progname);
      break;
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=argv[0];
  argc-=optind;
  argv+=optind;

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  inputfilename      =  *argv;
  parmtopfilename    =  *++argv;
  outputfilename1    =  *++argv;
  outputfilename2    =  *++argv;
  outputfilename3    =  *++argv;

  /************************************************************************/
  /* if ((inMODE==CA && outMODE==AA) || (inMODE==CA && outMODE==HV) 	  */
  /*     || (inMODE==HV && outMODE==AA) || (inMODE==HV && outMODE==CA) )  */
  /*   printf("error\n");						  */
  /************************************************************************/
  
  parmtop=efopen(parmtopfilename,"r");
  readParmtop(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  traj=mync_get_trj_aw(inputfilename,MODE,IOMODE,numatom,&numatomp,&numstep);
  //  traj=myncL_get_trj_aw_b(inputfilename,inMODE,outMODE,IOMODE,numatom,&numatomp,&numstep);
  cov        = (double *)gcemalloc(sizeof(double)*numatomp*3*numatomp*3);
  eigenvalue = (double *)gcemalloc(sizeof(double)*numatomp*3);

  pca_norm(traj,numstep,numatomp);
  pca_covm(traj,numstep,numatomp,cov);
  pca_diag(cov,eigenvalue,numatomp);
  pca_proj(traj,cov,numstep,numatomp);

  outputfile1=efopen(outputfilename1,"w");
  for (i=0;i<numstep;++i) {
    for (j=0;j<spdim;++j) {
      fprintf(outputfile1," %e",traj[i*numatomp*3+j]);
    }
    fprintf(outputfile1,"\n");
  }
  fclose(outputfile1);

  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<spdim;++i) {
    for (j=0;j<numatomp*3;++j) {
      fprintf(outputfile2," %e",cov[i*numatomp*3+j]);
    }
    fprintf(outputfile2,"\n");
  }
  fclose(outputfile2);

  outputfile3=efopen(outputfilename3,"w");
  io_outputdata_f(outputfile3,numatomp*3,eigenvalue);
  fclose(outputfile3);

  return 0;
}

void USAGE(char *progname) {
  printf("-b [CA]");
  printf("-c [exclude all H]");
  printf("-K [amber]");
  printf("%s inputfilename parmtopfilename outputfilename1(trj) outputfilename2(vec) outputfilename3(val)\n",progname);
}
