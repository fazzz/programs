
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "bestfit.h"
#include "bestfitL.h"
#include "PTL.h"
#include "EF.h"

#include "netcdf_mineL.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,s;
  int numatom,numstep,numite,numinterval=1;
  int MODE=AA,NCMODE=MD,INMODE=AA,flago='x';
  double epsilon=0.0001;

  int MAXC=20;

  double *mass,massCA,*coord_ref,*rmsd_trj;
  double rmsd_ave=0.0,rmsd_ave_old=0.0;
  char *inputfilename,*parmtopfilename,outputfilename[200];
  char *outputfilenamebase,outputfilenamermsd[200];
  char *progname;

  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;

  FILE *inputfile, *parmtop;
  FILE *outputfile,*outputfilermsd,*log;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  int opt_idx=1;

  struct option long_opt[] = {
    {"Amber",0,NULL,'A'},
    {"CA",0,NULL,'C'},
    {"H",0,NULL,'H'},
    {"INCA",0,NULL,'I'},
    {"INHV",0,NULL,'J'},
    {"epsilon",1,NULL,'e'},
    {"numite",1,NULL,'n'},
    {"interval",1,NULL,'i'},
    {"MAXC",1,NULL,'M'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"ACHIJhe:n:i:M:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'A':
      NCMODE=AA; break;
    case 'C':
      MODE=CA; break;
    case 'H':
      MODE=HV; break;
    case 'I':
      INMODE=CA; break;
    case 'J':
      INMODE=HV; break;
    case 'e':
      epsilon=atof(optarg); break;
    case 'n':
      numite=atoi(optarg);  break;
    case 'i':
      numinterval=atoi(optarg); break;
    case 'M':
      MAXC=atoi(optarg); break;
    case 'h':
      USAGE(progname);   break;
    default:
      USAGE(progname);  exit(1);
    }
  }

  progname=argv[0];
  argc-=optind;
  argv+=optind;

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  inputfilename      =  *argv;
  parmtopfilename    =  *++argv;
  outputfilenamebase =  *++argv;

  parmtop=efopen(parmtopfilename,"r");
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  if (INMODE==CA) {
    j=0;
    for (i=0;i<numatom;++i) {
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
	++j; massCA=AP.AMASS[i];
      }
    }
    numatom=j;
  }
  else if (INMODE==HV) {
    j=0;
    for (i=0;i<numatom;++i) {
      if (strncmp(AP.IGRAPH[i],"H",1)!=0){
	++j;
      }
    }
    numatom=j;
  }

  mass=(double *)gcemalloc(sizeof(double)*numatom);
  if (INMODE==CA) for (i=0;i<numatom;++i) mass[i] = massCA;
  else if (INMODE==HV) {
    j=0;
    for (i=0;i<AP.NATOM;++i) {
      if (strncmp(AP.IGRAPH[i],"H",1)!=0){ 
	mass[j] = AP.AMASS[i];
	++j;
      }
    }
  }
  else for (i=0;i<numatom;++i) mass[i] = AP.AMASS[i];

  if (NCMODE==AMBER) numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MD);
  else numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MCD);

  coord_ref=(double *)gcemalloc(sizeof(double)*numatom*3);
  rmsd_trj =(double *)gcemalloc(sizeof(double)*numstep);

  for (i=0;i<MAXC;++i) {
    if (i>0) sprintf(inputfilename,"%s_cyc_%d_bf.trj",outputfilenamebase,i);

    if (flago=='o') sprintf(outputfilename,"%s_bf.trj",outputfilenamebase);
    else sprintf(outputfilename,"%s_cyc_%d_bf.trj",outputfilenamebase,i+1);

    if (NCMODE==AMBER) numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MD);
    else numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);

    rmsd_trj =(double *)gcemalloc(sizeof(double)*numstep);

    if (i>0) numinterval=1;

    ave_coord_ncbL(coord_ref,inputfilename,numatom,numstep,numinterval,MODE,NCMODE);

    bf_trajectry_ncbL(numatom,numstep,numinterval,mass,coord_ref,rmsd_trj,inputfilename,outputfilename,MODE,flago,NCMODE,INMODE);
    if (flago=='o') {
      sprintf(outputfilenamermsd,"%s_rmsd.txt",outputfilenamebase);
      outputfilermsd=efopen(outputfilenamermsd,"w");
      fprintf(outputfilermsd,"step    rmsd\n");
      for (j=0;j<numstep;++j)
	fprintf(outputfilermsd,"%d %12.8lf\n",j+1,rmsd_trj[j]);
      fclose(outputfilermsd);
      break;
    }

    flago='x';
    rmsd_ave=0.0;
    for (j=0;j<numstep;++j) rmsd_ave=(j*rmsd_ave+rmsd_trj[j])/(j+1);
    printf("cyc: %d  rmsd: %10.6lf\n",i+1,rmsd_ave);
    if (abs(rmsd_ave-rmsd_ave_old) < epsilon && i > 0) flago='o';
    rmsd_ave_old=rmsd_ave;
  }
   return 0;
}

void USAGE(char *progname) {
  printf("-C [alpha carbon] \n");
  printf("-H [exclude hydrogen] \n");
  printf("-K [Amber type] \n");
  printf("-h [help] \n");
  printf("-e [epsilon criteria of delta rmsd] \n");
  printf("-n [num of iteration (<10)] \n");
  printf("%s inputfilename parmtopfilename outputfilenamebase\n",progname);
}
