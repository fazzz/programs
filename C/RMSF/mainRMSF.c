#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "RMSF.h"
#include "PT.h"
#include "IOV2.h"
#include "EF.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int inistep,lasstep,numstep,numatom;

  int MODE,netcdfflag;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double ***crd;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;

  double *rmsf;

  char *inputfilename,*outputfilename,*parmfilename;
  FILE *inputfile,*outputfile,*parmfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {"a",0,NULL,'a'},
    {"b",0,NULL,'b'},
    {"c",0,NULL,'c'},
    {"ini",1,NULL,'i'},
    {"las",1,NULL,'l'},
    {"netcdf",0,NULL,1},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"habci:l:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 1:
      netcdfflag=ON;
      break;
    case 'i':
      inistep=atoi(optarg);
      break;
    case 'n':
      lasstep=atoi(optarg);
      break;
    case 'a':
      MODE=AA;
      break;
    case 'b':
      MODE=CA;
      break;
    case 'c':
      MODE=HV;
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  parmfilename = *++argv;
  outputfilename = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
 
  if (MODE==AA) numatom=AP.NATOM;
  else if (MODE==CA) numatom=AP.NRES;
  else if (MODE==HV) numatom=PT_countatomtype("H",1,NO);

  crd=(double ***)gcemalloc(sizeof(double **)*numatom);
  for (i=0;i<3;++i)crd[i]=(double **)gcemalloc(sizeof(double *)*3);
  
  if (netcdfflag==OFF) {
    ;
  }
  else {
    inputfile=efopen(inputfilename,"r");
    numstep=IOV2_scantrj_wrange(inputfile,inistep,lasstep,crd,MODE);
    fclose(inputfile);
  }

  rmsf=RMSF_prot(crd,numatom,numstep);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numatom;++i) fprintf(outputfile,"%lf\n",rmsf[i]);
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-b] CA \n");
  printf("[-c] HV \n");
  printf("[-h] help \n");
  printf("[--netcdf] netcdf input \n");
  printf("%s [-b] [-c] [-h] inputfilename parmfilename outputfilename \n",progname);
}
