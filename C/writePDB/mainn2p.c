#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PDBL.h"

#include "TOPO.h"
#include "PT.h"
#include "EF.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0
#define SUM 1
#define SERIES 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0;
  int numstep,numatom,interval=1;

  int IHflag=OFF,AMBncflag=OFF,MODE=AA;
  int OUTMODE;

  PDBLF PDBL;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crd;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;

  char *inputfilename,*outputfilename,*parmfilename,*progname,*logfilename;
  char outputfilename2[10000];
  FILE *outputfile,*parmfile,*logfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"int",1,NULL,'i'},
    {"PDBs",1,NULL,'P'},
    {"CA",0,NULL,'C'},
    {"H",0,NULL,'H'},
    {"Amber",0,NULL,'K'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hCHKi:P:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'i':
      interval=atoi(optarg);
      break;
    case 'C':
      MODE=CA;
      break;
    case 'H':
      MODE=HV;
      break;
    case 'P':
      logfilename=optarg;
      OUTMODE=SERIES;
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
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  PDBL.numatom=AP.NATOM;
  PDBL.PDBLa=(PDBLA *)gcemalloc(sizeof(PDBLA)*AP.NATOM);
  readPDBLdatafromParmtop(PDBL);

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  if (AMBncflag==OFF) numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);
  else numstep=mync_get_present_step_AMBER(inputfilename,&nc_id_MD);

  if (OUTMODE!=SERIES) outputfile=efopen(outputfilename,"w");
  else if (OUTMODE==SERIES) logfile=efopen(logfilename,"w");

  for (i=0;i<numstep;++i) {
    if (AMBncflag==OFF) mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
    else  mync_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id_MD,crd_nc);

    if (i%interval == 0) {
      for (j=0;j<AP.NATOM;++j) for (k=0;k < 3; ++k) PDBL.PDBLa[j].coord[k]=crd_nc[j][k];

      if (OUTMODE!=SERIES) {
	fprintf(outputfile,"MODEL\n");
	//    writPDBL(outputfile,PDBL);
	writPDBL_wopt(outputfile,PDBL,MODE);
	fprintf(outputfile,"ENDMOD\n");
      }
      else if (OUTMODE==SERIES) {
	++l;
	sprintf(outputfilename2,"%s_%d.pdb",outputfilename,l);
	outputfile=efopen(outputfilename2,"w");
	writPDBL_wopt(outputfile,PDBL,MODE);
	fclose(outputfile);
	fprintf(logfile,"%s\n",outputfilename2);
      }
    }
  }
  if (OUTMODE!=SERIES) fclose(outputfile);
  fclose(logfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-H] (include Hydrogen) \n");
  printf("[--CA] (include Hydrogen) \n");
  printf("[--PDBs] logfilename \n");
  printf("[--Amber] (amber) \n");
  printf("[--int] interval \n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmfilename outputfilename \n",progname);
}
