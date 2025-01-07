
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "EF.h"
#include "IO.h"

#include "STRING.h"

#include "TOPO.h"
#include "RAND.h"
#include "BOXMULL.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0

#define kb 1.98723e-3

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  double *f;
  double **fb,**fa,**fd;
  double v,vext[2];
  double kd[2],n[2];
  double kapp;
  double dtheta[2];

  double pi;

  double eta,etaext;
  double beta=300.0,betaext=1000.0;
  double dt=0.01;
  int numstep=10000,intervat=1;
  double *crd,dihed1,dihed2;

  char *progname;
  char *inifilename,*outfilename,*log="log";
  FILE *inifile,*log;

  double crd_nc[5][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  pi=acos(-1.0);

  kd[0]=2.5;
  kd[1]=2.5;
  n[0]=3.0;
  n[1]=3.0;

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

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inifilename  = *argv;
  outfilename  = *++argv;

  beta=1.0/(kb*beta);
  betaext=1.0/(kb*betaext);

  crd=(double **)gcemalloc(sizeof(double *)*5);
  f=(double **)gcemalloc(sizeof(double *)*5);
  for (i=0;i<5;++i) {
    crd[i]=(double *)gcemalloc(sizeof(double)*3);
    f[i]=(double *)gcemalloc(sizeof(double)*3);
  }

  inifile=efopen(inifilename,"r");
  io_scanconf_Amber_ini(inifile,5,crd);
  fclose(inifile);

  for (i=0;i<numstep;++i) {
    z_string_FASYS_calcforce(crd,f,kd,n,fb,fa,fd);
    v=z_string_FASYS_calcpote(crd,kd,n);

    dihed1=dih(atom[0],atom[1],atom[2],atom[3]);
    dihed2=dih(atom[1],atom[2],atom[3],atom[4]);

    vext[0]=0.5*kapp*(z[0]-dihed1)*(z[0]-dihed1);
    vext[1]=0.5*kapp*(z[1]-dihed2)*(z[2]-dihed1);

    fext[0]=-kapp*(z[0]-dihed1);
    fext[1]=-kapp*(z[1]-dihed2);

    for (j=0;j<5*3;++j) {
      dihed1=dih(atom[0],atom[1],atom[2],atom[3]);
      dihed2=dih(atom[1],atom[2],atom[3],atom[4]);
      dtheta[0][j]=1.0/(2.0*h)*(dihed1ph-dihed1mh);
      dtheta[1][j]=1.0/(2.0*h)*(dihed1ph-dihed1mh);
    }

    for (j=0;j<5*3;++j) f[j]=fext[0]*dtheta[0][j]+fext[1]*dtheta[1][j];

    for (j=0;j<5*3;++j)
      crd[j]+=-dt/eta*f[j]+sqrt(2.0*dt/beta)/eta*Box_Muller(i,0.0,1.0);

    z[0]+=-dt/eta*fext[0]+sqrt(2.0*dt/betaext)/etaext*Box_Muller(i,0.0,1.0);
    z[1]+=-dt/eta*fext[1]+sqrt(2.0*dt/betaext)/etaext*Box_Muller(i,0.0,1.0);
    
    if (i%intervalt==0) {
      for (j=0;j<5;++j)
	for (k=0;k<3;++k) 
	  crd_nc[j][k]=crd[j*3+k];

      mync_put_crd_ene_MCD(nc_id_MCD,l,crd_nc,ene,drest);
      ++l;
    }
  }

  fclose(accplog);
  nc_close((nc_id_MCD.ncid));


  return 0;
}
 
void USAGE(char *progname) {
  printf("[-h] -- help\n");
  printf("USAGE: %s inifilename outfilename\n", progname);
}




