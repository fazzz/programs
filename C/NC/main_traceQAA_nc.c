
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "NC.h"
#include "PTL.h"
#include "EF.h"
#include "TOPO.h"
#include "MB.h"

#include "netcdf_mine.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numstep,HMODE=EXC;
  double QAA;

  double *crd,*crdref;
  int numatom,numres;

  double criteria=criteria_NC;
  int numncaa,numncres;
  int  *index_natatt,**ncmap;
  double atom1[3],atom2[3];

  double *cradii_natatt;
  double length;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_AMBER nc_id;

  char *progname;
  char *inputfilename,*parmtopname,*crdreffilename;
  char *outputfilename;

  FILE *crdfile,*parmtop,*crdreffile;
  FILE *outputfile;

  while((c=getopt(argc,argv,"hc:"))!=-1) {
    switch(c) {
    case 'c':
      criteria=atof(argv);
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
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdref=(double *)gcemalloc(sizeof(double)*numatom*3);

  crdreffile=efopen(crdreffilename,"r");
  io_scanconf_Amber_ini(crdreffile,numatom,crdref);
  fclose(crdreffile);

  numstep=mync_get_present_step_AMBER(inputfilename,&nc_id);

  ncmap=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) ncmap[i]=(int *)gcemalloc(sizeof(int)*numatom);
  index_natatt=make_native_contact_list_aa_2(&numncaa,&numncres,crdref,numatom,numres,criteria,ncmap,EXC);
  cradii_natatt=(double *)gcemalloc(sizeof(double)*numncaa);

  for (i=0;i<numncaa;++i) {
    for (j=0;j<3;++j) {
      atom1[j]=crdref[(index_natatt[i*2])*3+j];
      atom2[j]=crdref[(index_natatt[i*2+1])*3+j];
    }
    length=len(atom1,atom2);
    cradii_natatt[i]=len(atom1,atom2);
  }

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
      mync_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
      for (j=0;j<numatom;++j)for (k=0;k<3;++k)crd[j*3+k]=crd_nc[j][k];
      QAA=count_native_contact_aa(numncres,crd,numatom,numres,ncmap,cradii_natatt);

      fprintf(outputfile,"%e\n",QAA);
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmfilename outputfilename \n",progname);
}
