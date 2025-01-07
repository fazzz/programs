
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "EF.h"
#include "IO.h"

#include "SBFF.h"
#include "TOPO.h"
#include "RAND.h"
#include "BOXMULL.h"

#include "PT.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0

#define kb 1.98723e-3

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int d;
  int num_dec_temp=10;
  double f;
  int numarg;
  int numstep=10000,intervalt=1,intervale=1,intervald=1,num_dec_step;
  int numatom;

  double pi;

  double beta=300.0,deltaE;
  double dx=0.1;
  double *crd,*crd_trial,*crd1,*crd2;
  double *cv,*cv_trial,*cv_x,*cv_x_trial;
  double delta=1.0,deltaV;
  double p_t,p_t_trial;
  struct potential_SBAA ene1,ene1_trial,ene2,ene2_trial;

  char *progname;
  char *inifile1name,*inifile2name,*inifile3name,*parmtopname,*outfilename;
  FILE *inifile1,*inifile2,*inifile3,*parmtop;

  double /***crd_nc,*/*energy;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_SBAAMCD nc_id_MCD;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  pi=acos(-1.0);

  progname=argv[0];
  while((c=getopt(argc,argv,"chl:n:i:j:k:b:d:"))!=-1) {
    switch(c) {
    case 'h':
      USAGE(progname);
      exit(1);
    case 'n':
      numstep=atoi(optarg);
      break;
    case 'l':
      delta=atof(optarg);
      break;
    case 'i':
      intervalt=atoi(optarg);
      break;
    case 'j':
      intervale=atoi(optarg);
      break;
    case 'k':
      intervald=atoi(optarg);
      break;
    case 'b':
      beta=atof(optarg);
      break;
    case 'd':
      dx=atof(optarg);
      break;
    default:
      USAGE(progname);
      exit(1);
    }
  }

  argc-=optind;
  argv+=optind;

  numarg=4;

  if (argc < numarg) {
    USAGE(progname);
    exit(1);
  }
  inifile1name  = *argv;
  inifile2name  = *++argv;
  inifile3name  = *++argv;
  parmtopname  = *++argv;
  outfilename  = *++argv;

  beta=1.0/(kb*beta);

  parmtop=efopen(parmtopname,"r");
  readParmtop(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double )*numatom*3);
  crd_trial=(double *)gcemalloc(sizeof(double )*numatom*3);
  crd1=(double *)gcemalloc(sizeof(double )*numatom*3);
  crd2=(double *)gcemalloc(sizeof(double )*numatom*3);

  energy=(double *)gcemalloc(sizeof(double)*(NDIMS_ENERGY-1));

  inifile1=efopen(inifile1name,"r");
  io_scanconf_Amber_ini(inifile1,numatom,crd1);
  fclose(inifile1);

  inifile2=efopen(inifile2name,"r");
  io_scanconf_Amber_ini(inifile2,numatom,crd2);
  fclose(inifile2);

  inifile3=efopen(inifile3name,"r");
  io_scanconf_Amber_ini(inifile3,numatom,crd);
  fclose(inifile3);

  mync_create_def_SBAAMCD(outfilename,numatom,&nc_id_MCD);

  SBAAff_set_parameters_default(&ene1);
  SBAAff_set_parameters_default(&ene2);
  SBAAff_set_parameters_default(&ene1_trial);
  SBAAff_set_parameters_default(&ene2_trial);

  SBAAMBff_set_calcff(&ene1,&ene2,&deltaV,crd1,crd2,numatom);
  SBAAMBff_set_calcff(&ene1_trial,&ene2_trial,&deltaV,crd1,crd2,numatom);
  p_t=SBAAMBff_calcff(crd,numatom,&ene1,&ene2,delta,deltaV);
  p_t_trial=SBAAMBff_calcff(crd,numatom,&ene1_trial,&ene2_trial,delta,deltaV);

  for (i=0;i<numatom*3;++i)
    crd_trial[i]=crd[i];
  
  l=0;
  for (i=0;i<numstep;++i) {
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k)
	crd_trial[j*3+k]=crd[j*3+k]+dx*Box_Muller(i,0.0,1.0);
      p_t_trial=SBAAMBff_calcff(crd_trial,numatom,&ene1_trial,&ene2_trial,delta,deltaV);
      deltaE=p_t_trial-p_t;
      if((c=Metropolis(beta*deltaE))==1) {
	for (k=0;k<numatom*3;++k)
	  crd[k]=crd_trial[k];
	p_t=p_t_trial;
	memcpy(&ene1,&ene1_trial,sizeof(ene1_trial));
	memcpy(&ene2,&ene2_trial,sizeof(ene2_trial));
      }
    }

    if (i%intervalt==0) {
      for (j=0;j<numatom;++j)
	for (k=0;k<3;++k) 
	  crd_nc[j][k]=crd[j*3+k];

      mync_put_crd_ene_SBAAMCD(nc_id_MCD,l,crd_nc,p_t);
      ++l;
    }
  }
    
  nc_close((nc_id_MCD.ncid));

  return 0;
}
 
void USAGE(char *progname) {
  printf("[-h] -- help\n");
  printf("[-c] -- help for option\n");
  printf("[-a drestfilename ] -- drestflag\n");
  printf("[-n numstep ]\n");
  printf("[-i intervalt ]\n");
  printf("[-j intervale ]\n");
  printf("[-k intervald ]\n");
  printf("[-b beta ]\n");
  printf("[-d dx ]\n");
  printf("USAGE: %s crdfilename(ref1) crdfilename(ref2) inifilename parmtopname outfilename\n", progname);
}

