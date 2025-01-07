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

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int flagb=1,flaga=1,flagd=1;
  int numatom,numcys;
  int nbond,nangl,ndihed;

  double *crd,p,*f;
  double dx=0.001;

  int **atoms_b,**atoms_a,**atoms_d;
  double *K_bond,*K_angl,*V_dihe;
  double *eq_bond, *eq_angl,*n_dihe,*the_dihe;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*outputfilename,*parmfilename;
  FILE *inputfile,*outputfile,*parmfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {"b",0,NULL,'b'},
    {"a",0,NULL,'a'},
    {"d",0,NULL,'d'},
    {"x",0,NULL,'x'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hbadx:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'b':
      flagb=0;
      break;
    case 'a':
      flaga=0;
      break;
    case 'd':
      flagd=0;
      break;
    case 'x':
      dx=atof(optarg);
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

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  parmfilename = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);

  numcys=disulfid_count_numcys();

  atoms_b=(int **)gcemalloc(sizeof(int *)*numcys);
  atoms_a=(int **)gcemalloc(sizeof(int *)*numcys*/*3*/4);
  atoms_d=(int **)gcemalloc(sizeof(int *)*numcys*2);
  K_bond=(double *)gcemalloc(sizeof(double)*numcys);
  eq_bond=(double *)gcemalloc(sizeof(double)*numcys);
  K_angl=(double *)gcemalloc(sizeof(double)*numcys*4);
  eq_angl=(double *)gcemalloc(sizeof(double)*numcys*4);
  V_dihe=(double *)gcemalloc(sizeof(double)*numcys*2);
  n_dihe=(double *)gcemalloc(sizeof(double)*numcys*2);
  the_dihe=(double *)gcemalloc(sizeof(double)*numcys*2);

  for (i=0;i<numcys;++i) {
    atoms_b[i]=(int *)gcemalloc(sizeof(int)*3);
  }
  for (i=0;i<numcys*3;++i) {
    atoms_a[i]=(int *)gcemalloc(sizeof(int)*4);
  }
  for (i=0;i<numcys*2;++i) {
    atoms_d[i]=(int *)gcemalloc(sizeof(int)*5);
  }

  disulfid_read_parm(atoms_b,K_bond,eq_bond,&nbond,
		     atoms_a,K_angl,eq_angl,&nangl,
		     atoms_d,V_dihe,n_dihe,the_dihe,&ndihed);

  for (i=0;i<numcys;++i) {
    printf("K_b=%4.2lf eq_b=%4.2lf\n",K_bond[i],eq_bond[i]);
  }
  for (i=0;i<numcys*4;++i) {
    printf("K_a=%4.2lf eq_a=%4.2lf\n",K_angl[i],eq_angl[i]);
  }
  for (i=0;i<numcys*2;++i) {
    printf("V_d=%4.2lf n_d=%4.2lf the_d=%4.2lf\n",V_dihe[i],n_dihe[i],the_dihe[i]);
  }

  numatom=AP.NATOM;
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  f=(double *)gcemalloc(sizeof(double)*numatom*3);
  inputfile=efopen(inputfilename,"r");
  io_inputtrj_Amberform(inputfile,crd);
  fclose(inputfile);

  disulfid_calc_pf(crd,&p,f,
		   atoms_b,K_bond,eq_bond,nbond,flagb,
		   atoms_a,K_angl,eq_angl,nangl,flaga,
		   atoms_d,V_dihe,n_dihe,the_dihe,ndihed,flagd
		   );
  
  printf("p_ss=%4.2lf\n",p);  
  for (i=0;i<AP.NATOM;++i)
    if ( f[i*3+0]!=0.0 || f[i*3+0]!=0.0  || f[i*3+0]!=0.0  ) 
      printf("%d: %s %4.2lf %4.2lf %4.2lf\n",i+1,AP.IGRAPH[i],f[i*3+0],f[i*3+1],f[i*3+2]);

  disulfid_check_force_calc(crd,f,dx,
			    atoms_b,K_bond,eq_bond,nbond,flagb,
			    atoms_a,K_angl,eq_angl,nangl,flaga,
			    atoms_d,V_dihe,
			    n_dihe,the_dihe,ndihed,flagd);

  printf("\n");  
  for (i=0;i<AP.NATOM;++i)
    if ( f[i*3+0]!=0.0 || f[i*3+0]!=0.0  || f[i*3+0]!=0.0  ) 
      printf("%d: %s %4.2lf %4.2lf %4.2lf\n",i+1,AP.IGRAPH[i],f[i*3+0],f[i*3+1],f[i*3+2]);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmfilename outputfilename \n",progname);
}
