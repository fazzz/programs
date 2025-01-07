
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "EF.h"
#include "IO.h"
#include "TOPO.h"
#include "PT.h"
#include "netcdf_mine.h"

#define ON 1
#define OFF 0

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int numarg;

  int numatom;

  double *crd,atom_i[3],atom_j[3];
  double length;

  char *progname;
  char *inifilename,*parmtopname,*outfilename;
  FILE *inifile,*parmtop,*outfile;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
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

  argc-=optind;
  argv+=optind;

  numarg=3;

  if (argc < numarg) {
    USAGE(progname);
    exit(1);
  }
  inifilename  = *argv;
  parmtopname  = *++argv;
  outfilename  = *++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtop(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double )*numatom*3);

  inifile=efopen(inifilename,"r");
  io_scanconf_Amber_ini(inifile,numatom,crd);
  fclose(inifile);

  outfile=efopen(outfilename,"w");
  for (i=0;i<numatom;++i) {
    for (j=0;j<numatom;++j) {
      if (i<=j) {
	//	for (k=0;k<10;++k)
	  fprintf(outfile,"%5.3e ",0.0);
      }
      else {
	for (k=0;k<3;++k) {
	  atom_i[k]=crd[i*3+k];
	  atom_j[k]=crd[j*3+k];
	}

	length=len(atom_i,atom_j);
	if (length>3.0)
	  length=1.0;
	else
	  length=0.0;
	//	for (k=0;k<10;++k)
	  fprintf(outfile,"%5.3e ",length);
      }
      fprintf(outfile,"\n ");
    }
  }
  fclose(outfile);

}
 
void USAGE(char *progname) {
  printf("[-h] -- help\n");
  printf("USAGE: %s inifilename parmtopname outfilename\n", progname);
}

