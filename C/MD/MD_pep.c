
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "/home/yamamori/work/programs/yuMD2.1/src/massMD/src/disulfid.h"
//#include "/home/yamamori/work/programs/yuMD2.1/src/massMD/src/ParmTop.h"
#include "PTL.h"
#include "EF.h"

#include "FF.h"

#define ON 1
#define OFF 0

#define kb 1.98723e-3

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int d;
  double f;
  int drestflag=OFF;
  int numarg;
  int numstep=10000,intervalt=1,intervale=1,intervald=1;
  int numatom,numpara,numnb,num14;

  double pi;

  double beta=300.0,delta;
  double dx=0.1;
  double *crd,*crd_trial;
  double *cv,*cv_trial,*cv_x,*cv_x_trial;
  double drest,drest_trial;
  struct potential ene,ene_trial;
  struct drestparameters drestparm;

  char *progname;
  char *inifilename,*parmtopname,*indexfilename,*trjfilename,*enefilename,*enedfilename,*drestfilename;
  FILE *inifile,*parmtop,*indexfile,*trjfile,*enefile,*enedfile,*drestfile;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  int opt_idx=1;

  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  pi=acos(-1.0);

  progname=argv[0];
  while((c=getopt_long(argc,argv,"h",long_opt,&opt_idx))!=-1) {
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

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  inifilename  = *argv;
  parmtopname  = *++argv;
  enefilename  = *++argv;
  trjfilename  = *++argv;
  velfilename  = *++argv;

  beta=1.0/(kb*beta);

  parmtop=efopen(parmtopname,"r");
  readParmtop(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;
  numpara=AP.NTYPES*(AP.NTYPES+1)/2;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);

  inifile=efopen(inifilename,"r");
  io_scanconf_Amber_ini(inifile,numatom,crd);
  fclose(inifile);

  trjfile=efopen(trjfilename,"w");
  velfile=efopen(velfilename,"w");
  enefile=efopen(enefilename,"w");

  //  ff_set_calcffsp(&ene);

  //  ff_calcff(crd,numatom,&ene);

  disulfid_read_parm(atoms_b,K_bond,eq_bond,&nbond,
		     atoms_a,K_angl,eq_angl,&nangl,
		     atoms_d,V_dihe,n_dihe,the_dihe,&ndihed);

  numatom=AP.NATOM;
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  vel=(double *)gcemalloc(sizeof(double)*numatom*3);
  f=(double *)gcemalloc(sizeof(double)*numatom*3);
  inputfile=efopen(inputfilename,"r");
  io_inputtrj_Amberform(inputfile,crd);
  fclose(inputfile);

  for (i=0;i<numstep;++i) {
    disulfid_calc_pf(crd,&p,f,
		     atoms_b,K_bond,eq_bond,nbond,
		     atoms_a,K_angl,eq_angl,nangl,
		     atoms_d,V_dihe,n_dihe,the_dihe,ndihed
		     );
    for (i=0;i<numatom;++i) {
      vel
    
    }
    //      ff_calcff(crd_trial,numatom,&ene_trial);

    if (i%intervalt==0)
      io_outputconf(trjfile,numatom,crd,'x');
    if (i%intervald==0) {
      fprintf(enedfile,"%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",ene.p_t,ene.p_e_t,ene.p_LJ_t,ene.p_e_14_t,ene.p_LJ_14_t,ene.p_d_t,ene.p_a_t,ene.p_b_t);
    }
  }
    
  fclose(trjfile);
  fclose(velfile);
  fclose(enedfile);
  
  return 0;
}
 
void USAGE(char *progname) {
  printf("-h -- help\n");
  printf("USAGE: %s inifilename parmtopname enefilename trjfilename velfilename \n", progname);
}






