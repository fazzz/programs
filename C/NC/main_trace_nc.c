
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "NC.h"
#include "PT.h"
#include "EF.h"

#include "netcdf_mine.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numstep,HMODE=EXC;
  double q;

  int numnc,*indexncb;
  int **ncmap;
  double *crd,*crdref;
  int numatom,numres;
  double criteria=6.5;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  char *progname;
  char *inputfilename,*parmtopname,*crdreffilename;
  char *outputfilename;

  FILE *crdfile,*parmtop,*crdreffile;
  FILE *outputfile;

  while((c=getopt(argc,argv,"Hhc:"))!=-1) {
    switch(c) {
    case 'H':
      HMODE=INC;
      exit(1);
    case 'c':
      criteria=atof(optarg);
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
  inputfilename = *argv;
  parmtopname = *++argv;
  crdreffilename = *++argv;
  outputfilename = *++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtop(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdref=(double *)gcemalloc(sizeof(double)*numatom*3);

  crdreffile=efopen(crdreffilename,"r");
  io_scanconf_Amber_ini(crdreffile,numatom,crdref);
  fclose(crdreffile);

  numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);
  ncmap=(int **)gcemalloc(sizeof(int *)*numres);
  for (i=0;i<numres;++i) ncmap[i]=(int *)gcemalloc(sizeof(int)*numres);
  indexncb=make_native_contact_list(&numnc,crdref,numatom,numres,criteria,ncmap);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
      mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
      for (j=0;j<numatom;++j)
	for (k=0;k<3;++k) 
	  crd[j*3+k]=crd_nc[j][k];

      q=count_native_contact(numnc,crd,numatom,numres,indexncb,criteria,HMODE);

      fprintf(outputfile,"%e\n",q);
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmfilename outputfilename \n",progname);
}
