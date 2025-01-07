#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "ABA.h"
#include "GOLMAA.h"
#include "GOLMAA_set.h"

#include "PTL.h"
#include "EF.h"
#include "NC.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int numatom,numres,numstep;
  int interval=1;
  int **nb_matrix;
  double *frc,pot;

  double T=300;
  double k_B=1.98723e-3;

  double constant=1.0;
  
  int num_NC,**ncmap,**ncmap_aa;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_AMBER nc_id;

  double *crd,*refcrd,*mass;
  struct potential_GOLMAA e;
  double pot_d;
  double R_C_D=1.0;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*reffilename,*outputfilename,*parmfilename;
  FILE *inputfile,*reffile,*outputfile,*parmfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {"temp",1,NULL,'t'},
    {"int",1,NULL,'i'},
    {"rate",1,NULL,'|'},
    {"cons",1,NULL,'>'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"*t:i:|:>:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'i':
      interval=atoi(optarg);
      break;
    case 't':
      T=atof(optarg);
      break;
    case '|':
      R_C_D=atof(optarg);
      break;
    case '>':
      constant=atof(optarg);
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename     = *argv;
  reffilename       = *++argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];

  numstep=mync_get_present_step_AMBER(inputfilename,&nc_id);

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);
  reffile=efopen(reffilename,"r");
  getline(&line,&len,reffile);
  fscanf(reffile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile,"%lf",&refcrd[i*3+j]);
  fclose(reffile);

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  nb_matrix=(int **)gcemalloc(sizeof(int *)*numatom);
  GOLMAAff_set_calcff(&e,refcrd,numatom,nb_matrix,R_C_D,constant);

  outputfile=efopen(outputfilename,"w");

  for (i=0;i<numstep;++i) {
    mync_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
    for (j=0;j<numatom;++j)for (k=0;k<3;++k)crd[j*3+k]=crd_nc[j][k];
    
    GOLMAAff_calcff(crd,numatom,&e,OFF,ON,ON,nb_matrix);
    pot=e.p_t;

    fprintf(outputfile,"%10.8lf \n",pot/k_B/T);

  }
  fclose(outputfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename outputfilename\n",progname);
}


