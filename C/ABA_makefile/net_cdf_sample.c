
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "EF.h"

#define NXYZ 3
#define COORDINATE "coordinate"
#define REC_NAME "step"

#define UNITS "units"
#define TRJ_UNIT "angstroum"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int d;
  double f;
  int numarg;
  int numstep=10,intervalt=1;
  int numatom=2;

  double dx=0.1;
  double **crd;

  int ncid,xyz_dimid,molecule_dimid,rec_dimid;
  int trj_varid;
  int dimids[3];

  char *progname,*trjfilename="trj.nc";

  size_t start[3], count[3];

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  progname=argv[0];
  while((c=getopt(argc,argv,"chan:i:j:k:b:d:"))!=-1) {
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

  crd=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) crd[i]=(double *)gcemalloc(sizeof(double)*3);

  nc_create(trjfilename,NC_CLOBBER,&ncid);

  nc_def_dim(ncid,"xyz",NXYZ,&xyz_dimid);
  nc_def_dim(ncid,"molecule",numatom,&molecule_dimid);
  nc_def_dim(ncid,REC_NAME,NC_UNLIMITED,&rec_dimid);

  dimids[0] = rec_dimid;
  dimids[1] = molecule_dimid;
  dimids[2] = xyz_dimid;

  nc_def_var(ncid,COORDINATE,NC_DOUBLE,3,dimids,&trj_varid);
  nc_put_att_text(ncid,trj_varid,UNITS,strlen(TRJ_UNIT),TRJ_UNIT);

  nc_enddef(ncid);

  count[0] = 1;
  count[1] = numatom;
  count[2] = 3;
  start[1] = 0;
  start[2] = 0;

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      crd[i][j]=0.0;

  for (i=0;i<numstep;++i) {
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k)
	crd[j][k]=crd[j][k]+dx*Box_Muller(i,0.0,1.0);
    }

    if (i%intervalt==0) {
      start[0]=i;
      nc_put_vara_float(ncid,trj_varid,start,count,&crd);
    }
  }
    
  nc_close(ncid);
  
  return 0;
}
 
void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("USAGE: %s trjfilename \n", progname);
}

