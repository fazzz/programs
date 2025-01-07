#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "TOPO.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numstep,numatom;
  double length;

  int IHflag=OFF,crdflag=OFF;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double atom[2][3];
  double *crd;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  char *inputfilename,*outputfilename,*parmfilename,*progname;
  FILE *inputfile,*outputfile,*parmfile;

  while((c=getopt(argc,argv,"hHAs:"))!=-1) {
    switch(c) {
    case 's':
      numstep=atoi(optarg);
      break;
    case 'H':
      IHflag=ON;
      break;
    case 'A':
      crdflag=ON;
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

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  if (crdflag==OFF) numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);
  else {
    inputfile=efopen(inputfilename,"r");
    numstep=1;
  }
  outputfile=efopen(outputfilename,"w");

  for (i=0;i<numstep;++i) {
    if (crdflag==OFF) {
      mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd[j*3+k]=crd_nc[j][k];
    }
    else {
      d=io_scanconfwj(inputfile,numatom,crd,'x');
      if (d==-1)
	break;
      ++numstep;
    }
    if (IHflag==ON) {
      for (j=0;j<AP.NBONH;++j) {
	for (k=0;k<2;++k) for (l=0;l<3;++l) atom[k][l]=crd[(AP.BH[j][k])+l];
      
	length=len(atom[0],atom[1]);
	fprintf(outputfile,"%e ",length);
      }
    }

    for (j=0;j<AP.MBONA;++j) {
      for (k=0;k<2;++k) for (l=0;l<3;++l) atom[k][l]=crd[(AP.TA[j][k])+l];

	length=len(atom[0],atom[1]);
	fprintf(outputfile,"%e ",length);
    }
    fprintf(outputfile,"\n");
  }
  fclose(inputfile);
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-H] INCLUDE Hatom flag \n");
  printf("[-A] ON crdflag \n");
  printf("[-h] help \n");
  printf("%s [-H] [-h] inputfilename parmfilename outputfilename \n",progname);
}
