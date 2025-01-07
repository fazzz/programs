#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PDBL.h"

#include "PT.h"
#include "EF.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,n;
  int numatom,numPDB;
  int AMBncflag=OFF;

  PDBLF PDBL;
  
  extern char *optarg;
  extern int optind,opterr,optopt;

  int c;
  double *crd;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  struct potential e;

  char *inputfilename,*outputfilename,*parmfilename,*progname;
  char PDBfilename[1000];

  FILE *inputfile,*outputfile,*parmfile,*PDBfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"Amber",0,NULL,'K'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hK",long_opt,&opt_idx))!=-1) {
    switch(c) {
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

  myncL_create_def_MCD(outputfilename,numatom,&(nc_id_MCD));

  inputfile=efopen(inputfilename,"r");
  fscanf(inputfile,"%d",&numPDB);

  while ((c=getc(inputfile))!=-1){
    if (c!=' ' &&  c!='\n') {
      PDBfilename[0]=c;
      break;
    }
  }

  l=0;
  j=1;
  for (i=0;i<numPDB;++i) {
    while ((c=getc(inputfile))!=-1){
      if (c==' ' || c=='\n') {
	PDBfilename[j]='\0';
	PDBfile=efopen(PDBfilename,"r");
	readPDBL(PDBfile,PDBL,numatom);
	j=0;
	break;
      }
      else {
	PDBfilename[j]=c;
	++j;
      }
    }

    for (n=0;n<numatom;++n) for (k=0;k<3;++k) crd_nc[n][k]=PDBL.PDBLa[n].coord[k];
    myncL_put_crd_ene_MCD(nc_id_MCD,&l,crd_nc,e,0.0);
    ++l;
    
  }
  fclose(inputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[--Amber] (amber) \n");
  printf("[-h] help \n");
  printf("%s [--Amber] [-h] inputfilename parmfilename outputfilename \n",progname);
}
