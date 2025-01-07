#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "GOLMAA.h"
#include "GOLMAA_set.h"
#include "GOLMAA_check.h"
#include "GOLMAA_dbasin.h"

#include "PTL.h"
#include "EF.h"
#include "PDB.h"
#include "NC.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,ii,jj;
  int numatom,numres;
  int pdbflag=OFF;
  int flagd=ON,flagn=ON,flagr=ON;

  int flagfrccheck=OFF,flagadd=OFF;

  int **nb_matrix1,**nb_matrix2;

  int numspatom=11;
  double dx=0.00001;
  double *f,f_n[3],f_r[3];

  double delta=0.5,deltaV=0.5,ratio;
  
  double *crd,*refcrd1,*refcrd2;
  struct potential_GOLMAA_dbasin e;
  double R_C_D=1.0;

  PDBF PDB,PDBref;

  /****************************/
  /* int num_NCcrd,num_NCref; */
  /* int *ncmapcrd,*ncmapref; */
  /****************************/
  int numnc,*indexncb;
  int **ncmap;
  double Q_NC;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*reffilename1,*reffilename2,*parmfilename,*outputfilename;
  char *progname;
  FILE *inputfile,*parmfile,*reffile1,*reffile2,*outputfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"pdb",0,NULL,'p'},
    {"f",0,NULL,'f'},
    {"a",0,NULL,'a'},
    {"d",0,NULL,'d'},
    {"n",0,NULL,'n'},
    {"r",0,NULL,'r'},
    {"s",1,NULL,'s'},
    {"dx",1,NULL,'x'},
    {"rate",1,NULL,'|'},
    {"delta",1,NULL,'6'},
    {"deltaV",1,NULL,'V'},
    {"H",0,NULL,'H'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"pfadnrHs:|:x:6:V:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'p':
      pdbflag=ON;
      break;
    case 'f':
      flagfrccheck=ON;
      break;
    case 'a':
      flagadd=ON;
      break;
    case 'd':
      flagd=OFF;
      break;
    case 'x':
      dx=atof(optarg);
      break;
    case 's':
      numspatom=atoi(optarg);
      break;
    case '|':
      R_C_D=atof(optarg);
      break;
    case '6':
      delta=atof(optarg);
      break;
    case 'V':
      deltaV=atof(optarg);
      break;
    case 'H':
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  reffilename1 = *++argv;
  reffilename2 = *++argv;
  parmfilename = *++argv;
  outputfilename = *++argv;

  parmfile=efopen(parmfilename,"r");  
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;
  numres=AP.NRES;
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd1=(double *)gcemalloc(sizeof(double)*numatom*3);
  refcrd2=(double *)gcemalloc(sizeof(double)*numatom*3);
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

  reffile1=efopen(reffilename1,"r");
  if (pdbflag==OFF) {
    io_scanconf_Amber_ini(reffile1,numatom,refcrd1);
    fclose(reffile1);
  }
  else {
    readPDB(reffile1,PDB,numatom);
    for (i=0;i<numatom;++i)
      for (j=0;j<3;++j)
	refcrd1[i*3+j]=PDB.PDBa[i].coord[j];
  }

  reffile2=efopen(reffilename2,"r");
  if (pdbflag==OFF) {
    io_scanconf_Amber_ini(reffile2,numatom,refcrd2);
    fclose(reffile2);
  }
  else {
    readPDB(reffile2,PDB,numatom);
    for (i=0;i<numatom;++i)
      for (j=0;j<3;++j)
	refcrd2[i*3+j]=PDB.PDBa[i].coord[j];
  }

  GOLMAAff_dbasin_set_calcff(&e,refcrd1,refcrd2,numatom,R_C_D);

  if (flagadd==OFF)
    outputfile=efopen(outputfilename,"w");
  else
    outputfile=efopen(outputfilename,"a");

  ratio=GOLMAAff_dbasin_calcff_ratio(crd,numatom,&e,delta,deltaV);
  fprintf(outputfile,"p_t=%8.3lf ratio=%8.3lf\n",e.p_t,ratio);

  printf("p_t=%8.3lf ratio=%8.3lf\n",e.p_t,ratio);

  if (flagfrccheck==ON) {
    fprintf(outputfile,"f_x=%8.3lf f_y=%8.3lf f_z=%8.3lf \n",e.f_t[numspatom][0],e.f_t[numspatom][1],e.f_t[numspatom][2]);
  }

  if (flagfrccheck==ON) {
    f=GOLMAAff_dbasin_calcff_check(crd,numatom,numspatom,dx,refcrd1,refcrd2,R_C_D,delta,deltaV);
    printf("f_t_x=%8.3lf f_t_y=%8.3lf f_t_z=%8.3lf \n",e.f_t[numspatom][0],e.f_t[numspatom][1],e.f_t[numspatom][2]);
    printf("f_t_x=%8.3lf f_t_y=%8.3lf f_t_z=%8.3lf \n",f[0],f[1],f[2]);
  }

  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-help]\n");
  printf("%s inputfilename reffilename parmfilename\n",progname);
}

