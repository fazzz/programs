#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PTL.h"
#include "FFL.h"
#include "EF.h"
#include "RAND.h"
#include "BOXMULL.h"
#include "MD.h"

#include "netcdf_mineL.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,d;
  int numatom,numstep=10000,interval=100;
  double dt=0.001;
  double *frc,PE;

  int vMode=OFF;
  int bflag=ON,aflag=ON,dflag=ON,eflag=ON,LJflag=ON,e14flag=ON,LJ14flag=ON;

  double T0=300,T,K0,KE;
  double k_B=1.98723e-3;
  double UNITT=418.4070;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;

  double *crd,*mass,*vel,*vel0;

  struct potential e;
  struct force f;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*velfilename,*parmfilename;
  char *trjfilename,*outputfilename,*outputfilename2,*rstfilename="rstcrd",*rstvelfilename="rstvel";

  FILE *inputfilestat,*velfile,*parmfile;
  FILE *outputfile,*outputfile2,*rstfile,*rstvelfile;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"bnd",0,NULL,'b'},
    {"ang",0,NULL,'a'},
    {"dih",0,NULL,'d'},
    {"e14",0,NULL,'1'},
    {"l14",0,NULL,'4'},
    {"e",0,NULL,'e'},
    {"l",0,NULL,'l'},
    {"vMode",1,NULL,'v'},
    {"nums",1,NULL,'s'},
    {"temp",1,NULL,'t'},
    {"dt",1,NULL,'x'},
    {"int",1,NULL,'i'},
    {"rst",1,NULL,'{'},
    {"rstvel",1,NULL,'}'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"bad14elvhx:s:t:i:{:}:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'b':
      bflag=OFF;
      break;
    case 'a':
      aflag=OFF;
      break;
    case 'd':
      dflag=OFF;
      break;
    case '1':
      e14flag=OFF;
      break;
    case '4':
      LJ14flag=OFF;
      break;
    case 'e':
      eflag=OFF;
      break;
    case 'l':
      LJflag=OFF;
      break;
    case 's':
      numstep=atoi(optarg);
      break;
    case 't':
      T0=atof(optarg);
      break;
    case 'i':
      interval=atoi(optarg);
      break;
    case '{':
      rstfilename=optarg;
      break;
    case '}':
      rstvelfilename=optarg;
      break;
    case 'x':
      dt=atof(optarg);
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

  if (argc < 5) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  parmfilename      = *++argv;
  outputfilename    = *++argv;
  outputfilename2   = *++argv;
  trjfilename       = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i]=AP.AMASS[i];

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  vel=(double *)gcemalloc(sizeof(double)*numatom*3);
  vel0=(double *)gcemalloc(sizeof(double)*numatom*3);

  inputfilestat=efopen(inputfilename,"r");
  getline(&line,&len,inputfilestat);
  fscanf(inputfilestat,"%d",&d);
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(inputfilestat,"%lf",&crd[i*3+j]);
  fclose(inputfilestat);

  //  K0=0.5*(3*numatom-1.0)*k_B*T0;
  if ( vMode==OFF )
    K0=MD_Generate_inivelo(vel,mass,numatom,k_B*T0*UNITT);
  else {
    velfile=efopen(velfilename,"r");
    for (i=0;i<numatom;++i) for (j=0;j<3;++j) fscanf(velfile,"%lf",&vel[i*3+j]);
    fclose(velfile);
    K0=0.0;
    for (i=0;i<numatom;++i) for (j=0;j<3;++j) K0+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];
  }
  for (i=0;i>numatom*3;++i) vel0[i]=vel[i];
  T=K0/((3*numatom-1.0)*k_B)*2.0/UNITT;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  ffL_set_calcffandforce(&e,&f);
  ffL_calcffandforce(crd,numatom,&e,&f);

  myncL_create_def_MCD(trjfilename,numatom,&nc_id_MCD);
  outputfile=efopen(outputfilename,"w");
  outputfile2=efopen(outputfilename2,"w");
  for (i=0;i<numstep;++i) {
    KE=MD_Propagetor_Iso_JCP2003(crd,vel,mass,numatom,K0,dt,&e,&f);

    if (i%interval==0) {
      PE=e.p_b_t+e.p_a_t+e.p_d_t+0.5*e.p_e_t+0.5*e.p_LJ_t+0.5*e.p_e_14_t+0.5*e.p_LJ_14_t;
      KE=KE/UNITT;
      T=KE/((3*numatom-1.0)*k_B)*2.0;
	
      fprintf(outputfile,"%d %e %e %e %e\n",i+1,PE,KE,PE+KE,T);

      fprintf(outputfile2,"E_t    = %e \n",PE+KE);
      fprintf(outputfile2,"KE     = %e \n",KE);
      fprintf(outputfile2,"p_t    = %e \n",PE);
      fprintf(outputfile2,"p_es   = %e \n",e.p_e_t);
      fprintf(outputfile2,"p_LJ   = %e \n",e.p_LJ_t);
      fprintf(outputfile2,"p_14es = %e \n",e.p_e_14_t);
      fprintf(outputfile2,"p_14LJ = %e \n",e.p_LJ_14_t);
      fprintf(outputfile2,"p_dih  = %e \n",e.p_d_t);
      fprintf(outputfile2,"p_ang  = %e \n",e.p_a_t);
      fprintf(outputfile2,"p_bon  = %e \n",e.p_b_t);

      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crd[j*3+k];
      myncL_put_crd_ene_MCD(nc_id_MCD,l,crd_nc,e,0.0);
      ++l;
    }
  }
  fclose(outputfile);
  fclose(outputfile2);
  nc_close((nc_id_MCD.ncid));

  rstfile=efopen(rstfilename,"w");
  fprintf(rstfile,"ACE\n");
  fprintf(rstfile,"%d\n",&d);
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) fprintf(rstfile,"%e ",&crd[i*3+j]);
    fprintf(rstfile,"\n");
  }
  fclose(rstfile);

  rstvelfile=efopen(rstvelfilename,"w");
  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) fprintf(rstvelfile,"%e ",&vel[i*3+j]);
    fprintf(rstvelfile,"\n");
  }
  fclose(rstvelfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename parmfilename outputfilename\n",progname);
}


