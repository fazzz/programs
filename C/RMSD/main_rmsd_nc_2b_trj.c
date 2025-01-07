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
  int i,j,k,l,d;
  int initialstep,numstep,numatom,numatomp;
  double rmsd,*rot;
  double *mass;

  int AAflag=ON,CAflag=OFF,HVflag=OFF,crdflag=OFF,MODE,AMBERflag=OFF;
  
  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double *crd1,*crd2;
  double **crdA,**crdB;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD1,nc_id_MCD2;
  struct my_netcdf_out_id_AMBER nc_id_MD1,nc_id_MD2;

  char *inputfilename1,*inputfilename2,*outputfilename,*parmfilename;
  char *progname;
  FILE *outputfile,*parmfile,*condfile;

  while((c=getopt(argc,argv,"haAKbci:"))!=-1) {
    switch(c) {
    case 'i':
      initialstep=atoi(optarg);
      break;
    case 'a':
      MODE=AA;
      break;
    case 'A':
      crdflag=ON;
      break;
    case 'K':
      AMBERflag=ON;
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename1 = *argv;
  inputfilename2 = *++argv;
  parmfilename = *++argv;
  outputfilename = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  crd1=(double *)gcemalloc(sizeof(double)*numatom*3);
  crd2=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdA=(double **)gcemalloc(sizeof(double *)*3);
  for (j=0;j<3;++j) crdA[j]=(double *)gcemalloc(sizeof(double)*numatom);
  crdB=(double **)gcemalloc(sizeof(double *)*3);
  for (j=0;j<3;++j) crdB[j]=(double *)gcemalloc(sizeof(double)*numatom);
  rot=(double *)gcemalloc(sizeof(double)*9);
  mass=(double *)gcemalloc(sizeof(double)*numatom);

  if (crdflag==OFF) {
    if (AMBERflag==ON) {
      numstep=mync_get_present_step_AMBER(inputfilename1,&nc_id_MD1);
      numstep=mync_get_present_step_AMBER(inputfilename2,&nc_id_MD2);
    }
    else {
      numstep=mync_get_present_step_MCD(inputfilename1,&nc_id_MCD1);
      numstep=mync_get_present_step_MCD(inputfilename2,&nc_id_MCD2);
    }
  }

  outputfile=efopen(outputfilename,"w");

  for (i=0;i<numstep;++i) {
    if (AMBERflag==ON){
      mync_open_inq_get_sh_AMBER(inputfilename1,numatom,i,1,i+1,&nc_id_MD1,crd_nc);
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd1[j*3+k]=crd_nc[j][k];
      mync_open_inq_get_sh_AMBER(inputfilename2,numatom,i,1,i+1,&nc_id_MD2,crd_nc);
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd2[j*3+k]=crd_nc[j][k];
    }
    else {
      mync_open_inq_get_sh_MCD(inputfilename1,numatom,i,1,i+1,&nc_id_MCD1,crd_nc);
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd1[j*3+k]=crd_nc[j][k];
      mync_open_inq_get_sh_MCD(inputfilename2,numatom,i,1,i+1,&nc_id_MCD2,crd_nc);
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd2[j*3+k]=crd_nc[j][k];
    }

    l=0;
    for (j=0;j<numatom;++j) {
      if (MODE==CA) {
	if (strncmp(AP.IGRAPH[j],"CA",2)==0) {
	  mass[l]=AP.AMASS[j];
	  for (k=0;k<3;++k) {
	    crdA[k][l]=crd1[j*3+k];
	    crdB[k][l]=crd2[j*3+k];
	  }
	  ++l;
	}
      }
      else if (MODE==HV) {
	if (strncmp(AP.IGRAPH[j],"H",1)!=0) {
	  mass[l]=AP.AMASS[j];
	  for (k=0;k<3;++k) {
	    crdA[k][l]=crd1[j*3+k];
	    crdB[k][l]=crd2[j*3+k];
	  }
	  ++l;
	}
      }
      else {
	for (k=0;k<3;++k) {
	  crdA[k][l]=crd1[j*3+k];
	  crdB[k][l]=crd2[j*3+k];
	}
	++l;
      }
    }
    numatomp=l;

    fprintf(outputfile,"%lf\n",rmsd);    
    
    rmsd=CalcRMSDRotationalMatrix(crdB,crdA,numatomp,rot,mass);

  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-b] CA \n");
  printf("[-c] HV \n");
  printf("[-h] help \n");
  printf("%s [-b] [-c] [-h] inputfilename1 inputfilename2 parmfilename outputfilename \n",progname);
}
