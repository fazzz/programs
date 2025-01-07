#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "GOLM.h"
#include "GOLM_set.h"
#include "GOLM_check.h"

#include "PTL.h"
#include "EF.h"
#include "PDB.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int numatom;
  int pdbflag=OFF;
  int flagb=ON,flaga=ON,flagd=ON,flagn=ON,flagr=ON;

  int numspatom,numspres=2;
  double dx=0.000001;
  double *f;
  
  double *crd,*refcrd;
  struct potential_GOLM ene;

  PDBF PDB,PDBref;
  
  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*reffilename,*parmfilename;
  char *progname;
  FILE *inputfile,*parmfile,*reffile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"pdb",0,NULL,'p'},
    {"b",0,NULL,'b'},
    {"a",0,NULL,'a'},
    {"d",0,NULL,'d'},
    {"n",0,NULL,'n'},
    {"r",0,NULL,'r'},
    {"s",1,NULL,'s'},
    {"dx",1,NULL,'x'},
    {"H",0,NULL,'H'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"pbadnrHs:x:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'p':
      pdbflag=ON;
      break;
    case 'b':
      flagb=OFF;
      break;
    case 'a':
      flaga=OFF;
      break;
    case 'd':
      flagd=OFF;
      break;
    case 'n':
      flagn=OFF;
      break;
    case 'r':
      flagr=OFF;
      break;
    case 's':
      numspres=atoi(optarg);
      break;
    case 'x':
      dx=atof(optarg);
      break;
    case 'H':
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
  inputfilename = *argv;
  reffilename = *++argv;
  parmfilename = *++argv;

  parmfile=efopen(parmfilename,"r");  
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd=(double *)gcemalloc(sizeof(double)*numatom*3);
  if (pdbflag==ON) {
    PDB.PDBa=(PDBA *)gcemalloc(sizeof(PDBA)*numatom);
  }

  inputfile=efopen(inputfilename,"r");
  if (pdbflag==OFF) {
    io_scanconf_Amber_ini(inputfile,numatom,crd);
    fclose(inputfile);
  }
  else {
    readPDB(inputfile,PDB,numatom);
    for (i=0;i<numatom;++i)
      for (j=0;j<3;++j)
	crd[i*3+j]=PDB.PDBa[i].coord[j];
  }

  reffile=efopen(reffilename,"r");
  if (pdbflag==OFF) {
    io_scanconf_Amber_ini(reffile,numatom,refcrd);
    fclose(reffile);
  }
  else {
    readPDB(reffile,PDB,numatom);
    for (i=0;i<numatom;++i)
      for (j=0;j<3;++j)
	refcrd[i*3+j]=PDB.PDBa[i].coord[j];
  }

  GOLMff_set_calcff(&ene,refcrd,numatom);

  numspatom=PTL_res_ca(numspres);
  if (numspatom==-1) {
    printf("error of res num\n");
    exit(1);
  }

  GOLMff_calcff(crd,numatom,&ene,flagb,flaga,flagd,flagn,flagr);
  printf("p_t=%8.3lf \n",ene.p_t);
  printf("p_b=%8.3lf p_a=%8.3lf  p_d=%8.3lf \n",ene.p_b_t,ene.p_a_t,ene.p_d_t);
  printf("p_n=%8.3lf p_r=%8.3lf \n",ene.p_natatt_t,ene.p_repul_t);
  printf("f_x=%8.3lf f_y=%8.3lf f_z=%8.3lf \n",ene.f_t[numspatom][0],ene.f_t[numspatom][1],ene.f_t[numspatom][2]);


  f=GOLMff_calcff_check(crd,numatom,&ene,flagb,flaga,flagd,flagn,flagr,numspatom,dx);

  printf("f_x=%8.3lf f_y=%8.3lf f_z=%8.3lf \n",f[0],f[1],f[2]);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-help]\n");
  printf("%s inputfilename reffilename parmfilename\n",progname);
}

