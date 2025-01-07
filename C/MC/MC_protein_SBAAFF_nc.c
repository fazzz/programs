
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

  double accp,accn=0.0;

  double beta=300.0,deltaE;
  double dx=0.1;
  double *crd,*crd_trial,*crd_ref;
  double p_t,p_t_trial;
  struct potential_SBAA ene,ene_trial;

  char *progname;
  char *inifilename,*reffilename,*parmtopname,*outfilename,*accplogname="accp.log";
  FILE *inifile,*reffile,*parmtop,*accplog;

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
  while((c=getopt(argc,argv,"chn:i:j:k:b:d:a:"))!=-1) {
    switch(c) {
    case 'h':
      USAGE(progname);
      exit(1);
    case 'n':
      numstep=atoi(optarg);
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
    case 'a':
      accplogname=optarg;
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
  reffilename  = *argv;
  inifilename  = *++argv;
  parmtopname  = *++argv;
  outfilename  = *++argv;

  beta=1.0/(kb*beta);

  parmtop=efopen(parmtopname,"r");
  readParmtop(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double )*numatom*3);
  crd_trial=(double *)gcemalloc(sizeof(double )*numatom*3);
  crd_ref=(double *)gcemalloc(sizeof(double )*numatom*3);

  energy=(double *)gcemalloc(sizeof(double)*(NDIMS_ENERGY-1));

  reffile=efopen(reffilename,"r");
  io_scanconf_Amber_ini(reffile,numatom,crd_ref);
  fclose(reffile);

  inifile=efopen(inifilename,"r");
  io_scanconf_Amber_ini(inifile,numatom,crd);
  fclose(inifile);

  mync_create_def_SBAAMCD(outfilename,numatom,&nc_id_MCD);

  SBAAff_set_parameters_protein_default(&ene);
  SBAAff_set_parameters_protein_default(&ene_trial);

  SBAAff_set_protein_calcff(&ene,crd,numatom);
  SBAAff_set_protein_calcff(&ene_trial,crd,numatom);

  p_t=SBAAff_calcff(crd,numatom,&ene);
  p_t_trial=SBAAff_calcff(crd,numatom,&ene_trial);

  for (i=0;i<numatom*3;++i) crd_trial[i]=crd[i];
  
  accplog=efopen(accplogname,"w");
  fprintf(accplog,"#        a        p/n");
  l=0;
  for (i=0;i<numstep;++i) {
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k)
	crd_trial[j*3+k]=crd[j*3+k]+dx*Box_Muller(i,0.0,1.0);
    }
      SBAAff_calcff(crd_trial,numatom,&ene_trial);
      p_t_trial=ene_trial.p_t;
      deltaE=p_t_trial-p_t;
      if((c=Metropolis(beta*deltaE))==1) {
	accn=accn+c;
	for (k=0;k<numatom*3;++k) crd[k]=crd_trial[k];
	p_t=p_t_trial;
	memcpy(&ene,&ene_trial,sizeof(ene_trial));
      }
      //    }

    if (i%intervalt==0) {
      for (j=0;j<numatom;++j)
	for (k=0;k<3;++k) 
	  crd_nc[j][k]=crd[j*3+k];

      accp=accn/(i+1);
      fprintf(accplog,"%8.3d %8.3d %8.3lf\n",i+1,c,accp);
      mync_put_crd_ene_SBAAMCD(nc_id_MCD,l,crd_nc,p_t);
      ++l;
    }
  }
    
  fclose(accplog);
  nc_close((nc_id_MCD.ncid));

  return 0;
}
 
void USAGE(char *progname) {
  printf("[-h] -- help\n");
  printf("[-a drestfilename ] -- drestflag\n");
  printf("[-n numstep ]\n");
  printf("[-i intervalc ]\n");
  printf("[-j intervale ]\n");
  printf("[-k intervald ]\n");
  printf("[-b beta ]\n");
  printf("[-a accplogname ]\n");
  printf("[-d dx ]\n");
  printf("USAGE: %s crdfilename(ref) inifilename parmtopname outfilename\n", progname);
}


