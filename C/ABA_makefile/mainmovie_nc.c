#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "PDBL.h"

#include "TOPO.h"
#include "PT.h"
#include "EF.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int numstep,numatom;

  int IHflag=OFF,AMBncflag=OFF,MODE=AA;

  PDBLF PDBL;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crd;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;

  char *inputfilename,*outputfilename,*parmfilename,*progname;
  FILE *outputfile,*parmfile;

  while((c=getopt(argc,argv,"hCHKs:"))!=-1) {
    switch(c) {
    case 's':
      numstep=atoi(optarg);
      break;
    case 'C':
      //      IHflag=ON;
      MODE=CA;
      break;
    case 'H':
      //      IHflag=ON;
      MODE=HV;
      break;
    case 'K':
      AMBncflag=ON;
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
  numatom=AP.NATOM;

  PDBL.numatom=AP.NATOM;
  PDBL.PDBLa=(PDBLA *)gcemalloc(sizeof(PDBLA)*AP.NATOM);
  readPDBLdatafromParmtop(PDBL);

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  if (AMBncflag==OFF)
    numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);
  else
    numstep=mync_get_present_step_AMBER(inputfilename,&nc_id_MD);
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    if (AMBncflag==OFF)
      mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
    else
      mync_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id_MD,crd_nc);

    for (j=0;j<AP.NATOM;++j) for (k=0;k < 3; ++k) PDBL.PDBLa[j].coord[k]=crd_nc[j][k];

    fprintf(outputfile,"MODEL\n");
    //    writPDBL(outputfile,PDBL);
    writPDBL_wopt(outputfile,PDBL,MODE);
    fprintf(outputfile,"ENDMOD\n");
    
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-H] include Hydrogen \n");
  printf("[-K] amber \n");
  printf("[-h] help \n");
  printf("%s [-H] [-K] [-h] inputfilename parmfilename outputfilename \n",progname);
}
