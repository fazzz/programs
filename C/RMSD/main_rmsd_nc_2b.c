#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>


#include "rmsd.h"
#include "PT.h"
#include "EF.h"

#include "netcdf_mine.h"


#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int numatom,numatomp;
  double rmsd,*rot;
  double *mass;

  int AAflag=ON,CAflag=OFF,HVflag=OFF,MODE;
  
  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crd,*crdref;
  double **crdA,**crdB;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  char *inputfilename,*outputfilename,*parmfilename,*reffilename;
  char *progname;
  FILE *inputfile,*outputfile,*parmfile,*reffile;

  while((c=getopt(argc,argv,"habc"))!=-1) {
    switch(c) {
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
  reffilename = *++argv;
  parmfilename = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdref=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdA=(double **)gcemalloc(sizeof(double *)*3);
  for (j=0;j<3;++j) crdA[j]=(double *)gcemalloc(sizeof(double)*numatom);
  crdB=(double **)gcemalloc(sizeof(double *)*3);
  for (j=0;j<3;++j) crdB[j]=(double *)gcemalloc(sizeof(double)*numatom);
  rot=(double *)gcemalloc(sizeof(double)*9);
  mass=(double *)gcemalloc(sizeof(double)*numatom);

  inputfile=efopen(inputfilename,"r");
  io_scanconf_Amber_ini(inputfile,numatom,crd);
  fclose(inputfile);

  reffile=efopen(reffilename,"r");
  io_scanconf_Amber_ini(reffile,numatom,crdref);
  fclose(reffile);

  l=0;
  for (j=0;j<numatom;++j) {
    if (MODE==CA) {
      if (strncmp(AP.IGRAPH[j],"CA",2)==0) {
	mass[l]=AP.AMASS[j];
	for (k=0;k<3;++k) {
	  crdA[k][l]=crd[j*3+k];
	  crdB[k][l]=crdref[j*3+k];
	}
	++l;
      }
    }
    else if (MODE==HV) {
      if (strncmp(AP.IGRAPH[j],"H",1)!=0) {
	mass[l]=AP.AMASS[j];
	for (k=0;k<3;++k) {
	  crdA[k][l]=crd[j*3+k];
	  crdB[k][l]=crdref[j*3+k];
	}
	++l;
      }
    }
    else {
      for (k=0;k<3;++k) {
	crdA[k][l]=crd[j*3+k];
	crdB[k][l]=crdref[j*3+k];
      }
      ++l;
    }
  }
  numatomp=l;
  
  rmsd=CalcRMSDRotationalMatrix(crdB,crdA,numatomp,rot,mass);
  printf("%10.8e\n",rmsd);
  //  rmsd=rmsd_qcp(crd,crdref,numatom,MODE);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-b] CA \n");
  printf("[-c] HV \n");
  printf("[-h] help \n");
  printf("%s [-b] [-c] [-h] inputfilename reffilename parmfilename\n",progname);
}
