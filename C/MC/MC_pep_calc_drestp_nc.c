
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "PT.h"

#include "TOPO.h"
#include "EF.h"
#include "IO.h"

#include "netcdf_mine.h"

struct drestparameters {
  int numdrest;
  double k;
  double *dihed_equ;

  int *dihed_index;
};

void USAGE(char *progname);
double calc_direst(double *crd,int numatom,struct drestparameters drestparm, double *cv, double *cv_x);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int d;
  double f;
  int numarg;
  int numstep;
  int numatom;

  double pi;

  double *cv,*cv_x;
  double drest;
  struct drestparameters drestparm;

  char *progname;
  char *trjfilename,*parmtopname,*enefilename,*drestfilename;
  FILE *parmtop,*enefile,*drestfile;

  double crd[MAXATOM][3],*crd_data;
  struct my_netcdf_out_id_MCD nc_id_MCD;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  pi=acos(-1.0);

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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  enefilename  = *argv;
  trjfilename  = *++argv;
  drestfilename = *++argv;

  numatom=mync_get_numatom_MCD(trjfilename,&nc_id_MCD);
  crd_data=(double *)gcemalloc(sizeof(double)*numatom*3);
  numstep=mync_get_present_step_MCD(trjfilename,&nc_id_MCD);

  drestfile=efopen(drestfilename,"r");
  fscanf(drestfile,"%lf",&drestparm.k);
  fscanf(drestfile,"%d",&drestparm.numdrest);
  cv=(double *)gcemalloc(sizeof(double)*drestparm.numdrest);
  cv_x=(double *)gcemalloc(sizeof(double)*drestparm.numdrest);
  drestparm.dihed_index=(int *)gcemalloc(sizeof(int)*drestparm.numdrest*4);
  for (i=0;i<drestparm.numdrest;++i) {
    for (j=0;j<4;++j) {
      fscanf(drestfile,"%d",&d);
      drestparm.dihed_index[i*4+j]=d-1;
    }
    fscanf(drestfile,"%lf",&f);
    cv[i]=f/180*pi;
    if (cv[i]<-pi)
      cv[i]+=2.0*pi;
    else if (cv[i]>pi)
      cv[i]-=2.0*pi;
  }
  fclose(drestfile);

  enefile=efopen(enefilename,"w");
  fprintf(enefile,"#  Erest \n");
  for (i=0;i<numstep;++i) {
    mync_open_inq_get_sh_MCD(trjfilename,numatom,i,1,i+1,&nc_id_MCD,crd);
    l=0;
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
	crd_data[l]=crd[j][k];
	++l;
      }
    }
    drest=calc_direst(crd,numatom,drestparm,cv,cv_x);
    fprintf(enefile,"%d %12.8lf\n",i+1,drest);
  }
  fclose(enefile);

  return 0;
}

double calc_direst(double *crd,int numatom,struct drestparameters drestparm, double *cv, double *cv_x){
  int i,j,k;
  int numdrest;
  double atom[4][3],dihedang,dang;
  double drest=0.0,pi,f;

  pi=acos(-1.0);
  numdrest=drestparm.numdrest;

  for (i=0;i<numdrest;++i) {
    for (j=0;j<4;++j)
      for (k=0;k<3;++k)
	atom[j][k]=crd[drestparm.dihed_index[i*4+j]*3+k];
    /*cv_x[i]*/f = dih(atom[0],atom[1],atom[2],atom[3]);
    cv_x[i]=f;

    if ((dang=cv_x[i]-cv[i])>pi)
      dang-=2.0*pi;
    else if (dang<-1.0*pi)
      dang+=2.0*pi;

    drest+=0.5*(drestparm.k)*dang*dang;
  }

  return drest;
}

void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("-n numstep\n");
  printf("USAGE: %s parmtopname enefilename trjfilename drestfilename \n", progname);
}
