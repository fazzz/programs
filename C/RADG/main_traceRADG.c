
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "RADG.h"

#include "PTL.h"
#include "netcdf_mine.h"

#define ON 1
#define OFF 0

#define AA 0
#define CA 1
#define HA 2

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d;
  int numstep;
  double RADG,aRADG=0.0,vRADG=0.0;

  double **crd,*mass;
  int numatom;

  int AMBERMODEflag=OFF;
  int MODE=HA;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id;

  char *line;
  size_t le=0;

  char *progname;
  char *inputfilename,*parmtopname,*crdreffilename;
  char *outputfilename;

  FILE *crdfile,*parmtop,*crdreffile;
  FILE *outputfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"Amber",0,NULL,'A'},
    {"CA",0,NULL,'1'},
    {"HA",0,NULL,'2'},
    {"AA",0,NULL,'3'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"h123A:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'A':
      AMBERMODEflag=ON;
      break;
    case '1':
      MODE=CA;
      break;
    case '2':
      MODE=HA;
      break;
    case '3':
      MODE=AA;
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

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  parmtopname = *++argv;
  outputfilename = *++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  if (MODE==CA) {
    j=0;
    for (i=0;i<numatom;++i) {
      if (strncmp(AP.IGRAPH[i],"CA",2) == 0) {
	++j;
      }
    }
    numatom=j;
  }
  else if (MODE==HA) {
    j=0;
    for (i=0;i<numatom;++i) {
      if (strncmp(AP.IGRAPH[i],"H",1) != 0) {
	++j;
      }
    }
    numatom=j;
  }

  if (MODE==CA) {
    crd=(double **)gcemalloc(sizeof(double *)*numatom);
    for (i=0;i<numatom;++i) {
      crd[i]=(double *)gcemalloc(sizeof(double)*3);
    }
    mass=(double *)gcemalloc(sizeof(double)*numatom);
    j=0;
    for (i=0;i<AP.NATOM;++i) {
      if (strncmp(AP.IGRAPH[i],"CA",2) == 0) {
	mass[j]=AP.AMASS[i];
	++j;
      }
    }
  }
  else if (MODE==HA) {
    crd=(double **)gcemalloc(sizeof(double *)*numatom);
    for (i=0;i<numatom;++i) {
      crd[i]=(double *)gcemalloc(sizeof(double)*3);
    }
    mass=(double *)gcemalloc(sizeof(double)*numatom);
    j=0;
    for (i=0;i<AP.NATOM;++i) {
      if (strncmp(AP.IGRAPH[i],"H",1) != 0) {
	mass[j]=AP.AMASS[i];
	++j;
      }
    }
  }
  else {
    crd=(double **)gcemalloc(sizeof(double *)*numatom);
    for (i=0;i<numatom;++i) {
      crd[i]=(double *)gcemalloc(sizeof(double)*3);
    }
    mass=(double *)gcemalloc(sizeof(double)*numatom);
    for (i=0;i<numatom;++i) {
      mass[i]=AP.AMASS[i];
    }
  }

  if (AMBERMODEflag==ON) numstep=mync_get_present_step_AMBER(inputfilename,&nc_id);
  else numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    if (AMBERMODEflag==ON) mync_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id,crd_nc);
    else mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
    for (j=0;j<numatom;++j)for (k=0;k<3;++k)crd[j][k]=crd_nc[j][k];

    RADG=RADG_calc_radg(crd,mass,numatom);
    
    aRADG=(i*aRADG+RADG)/(i+1);
    vRADG=(i*vRADG+RADG*RADG)/(i+1);
    fprintf(outputfile,"%e\n",RADG);
  }
  fclose(outputfile);

  vRADG=sqrt(vRADG-aRADG*aRADG);

  printf("%s \n",inputfilename);
  printf("%s \n",outputfilename);
  printf("ave= %e var= %e\n",aRADG,vRADG);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("[--Amber] \n");
  printf("[--CA] \n");
  printf("[--HA] \n");
  printf("[--AA] \n");
  printf("%s [-h] inputfilename parmfilename outputfilename \n",progname);
}
