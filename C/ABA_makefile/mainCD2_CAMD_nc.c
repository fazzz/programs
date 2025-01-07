#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "netcdf_mine.h"

#include "PDB.h"
#include "MB.h"
#include "PTL.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

int usage(void);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numatom,numCAatom,numstep;
  int interval=1;

  double atom[4][3];
  double theta;
  double pi;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3],*crd;
  struct my_netcdf_out_id_MCD nc_id;

  char *progname;
  char *inputfilename,*outputfilename,*parmfilename;
  FILE *inputfile,*outputfile,*parmfile;

  while((c=getopt(argc,argv,"h"))!=-1) {
    switch(c) {
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
  parmfilename  =  *++argv;
  outputfilename = *++argv;

  pi=acos(-1.0);

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  j=0;
  for (i=0;i<numatom;++i) if (strncmp(AP.IGRAPH[i],"CA",2)==0) ++j;
  numCAatom=j;

  numstep=mync_get_present_step_MCD(inputfilename,&nc_id);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    mync_open_inq_get_sh_MCD(inputfilename,numCAatom,i,1,i+1,&nc_id,crd_nc);

    if ( (i%interval) == 0 ) {
      for (j=0;j<numCAatom-3;++j) {
	for (k=0;k<3;++k) {
	  atom[0][k]=crd_nc[j][k];
	  atom[1][k]=crd_nc[j+1][k];
	  atom[2][k]=crd_nc[j+2][k];
	  atom[3][k]=crd_nc[j+3][k];
	}
	
	theta=pick_dihed(atom[0],atom[1],atom[2],atom[3],0,0);
	if (theta > pi)  theta = theta -2.0*pi;
	else if (theta < -pi)  theta = theta+2.0*pi;
	fprintf(outputfile,"%e ",theta*180/acos(-1.0));
      }
      fprintf(outputfile,"\n");
    }
  }
  fclose(outputfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("%s inputfilename parmfilename outputfilename \n",progname);
}

 
