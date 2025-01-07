#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "GOLMAA_PROTEINS2008_set.h"

#include "PTL.h"
#include "EF.h"
#include "NC.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int resi,resj;
  int numatom,numres;
  int numncaa,numncres;

  int **ncmap,**ncmapres;
  double *crd,*refcrd;

  int nibnum=3,criteria=criteria_PROTEINS2008;

  struct potential_GOLMAA_PROTEINS2008 e_GOLM;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *reffilename,*parmfilename,*mapfilename;
  FILE *reffile,*parmfile,*mapfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"cutoff",1,NULL,'c'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hc:",long_opt,&opt_idx))!=-1) {
    switch(c) {
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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  reffilename       = *argv;
  parmfilename      = *++argv;
  mapfilename       = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);

  reffile=efopen(reffilename,"r");
  getline(&line,&len,reffile);
  fscanf(reffile,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(reffile,"%lf",&refcrd[i*3+j]);
  fclose(reffile);

  ncmap=GOLMAA_PROTEINS2008_ff_set_make_native_contact(refcrd,criteria,&numncaa,numatom,numres,nibnum);
        
  ncmapres=(int **)gcemalloc(sizeof(int *)*numres);
  for (i=0;i<numres;++i) ncmapres[i]=(int *)gcemalloc(sizeof(int)*numres);
  for (i=0;i<numres;++i) for (j=0;j<numres;++j) ncmapres[i][j]=0;

  numncres=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap[i][j]==0 ) {
	resi=PTL_resnum(i)-1;
	resj=PTL_resnum(j)-1;
	ncmapres[resi][resj]+=1;
	if (ncmapres[resi][resj]==1) ++numncres;
      }
    }
  }

  mapfile=efopen(mapfilename,"w");
  for (i=0;i<numres;++i)
    for (j=0;j<numres;++j)
      fprintf(mapfile,"%2d ",ncmapres[i][j]);
  fclose(mapfile);

  printf("%d %d",numncaa,numncres);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] refcrdfilename parmfilename mapfilename\n",progname);
}


