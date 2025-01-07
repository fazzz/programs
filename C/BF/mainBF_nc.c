
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "QUA.h"
#include "bestfit.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"
#include <netcdf.h>

#include "f2c.h"
#include "clapack.h"

#include "netcdf_mine.h"

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,s;
  int numatom,numstep,numite,numinterval=1;
  int MODE=AA,IOMODE=MD,flago='x';
  double epsilon=0.0001;

  double *mass,*coord_ref,*rmsd_trj;
  double rmsd_ave=0.0,rmsd_ave_old=0.0;
  char *inputfilename,*parmtopfilename,outputfilename[/*100*/200];
  char *outputfilenamebase,outputfilenamermsd[/*100*/200];
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

  while((c=getopt(argc,argv,"ACHKhe:n:i:"))!=-1) {
    switch(c) {
    case 'A':
      MODE=AA;
      break;
    case 'C':
      MODE=CA;
      break;
    case 'H':
      MODE=HV;
      break;
    case 'K':
      IOMODE=AMBER;
      break;
    case 'e':
      epsilon=atof(optarg);
      break;
    case 'n':
      numite=atoi(optarg);
      break;
    case 'i':
      numinterval=atoi(optarg);
      break;
    case 'h':
      USAGE(progname);
      break;
    default:
      USAGE(progname);
      exit(1);
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
  readParmtop(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  mass=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) mass[i] = AP.AMASS[i];

  if (IOMODE==AMBER)
    numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MD);
  else 
    numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);

  coord_ref=(double *)gcemalloc(sizeof(double)*numatom*3);
  rmsd_trj =(double *)gcemalloc(sizeof(double)*numstep);

  for (i=0;i<20/*MAXNUNITERATION*/;++i) {
    if (i>0) sprintf(inputfilename,"%s_cyc_%d_bf.trj",outputfilenamebase,i);

    if (flago=='o') {
      sprintf(outputfilename,"%s_bf.trj",outputfilenamebase);
    }
    else {
      sprintf(outputfilename,"%s_cyc_%d_bf.trj",outputfilenamebase,i+1);
    }

    if (IOMODE==AMBER)
      numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MD);
    else 
      numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);

    rmsd_trj =(double *)gcemalloc(sizeof(double)*numstep);

    if (i>0) numinterval=1;

    ave_coord_nc(coord_ref,inputfilename,numatom,numstep,MODE,IOMODE);
    //    ave_coord_ncb(coord_ref,inputfilename,numatom,numstep,numinterval,MODE,IOMODE);

    //    bf_trajectry_nc(numatom,numstep,mass,coord_ref,rmsd_trj,inputfilename,outputfilename,MODE,flago,IOMODE);
    bf_trajectry_ncb(numatom,numstep,numinterval,mass,coord_ref,rmsd_trj,inputfilename,outputfilename,MODE,flago,IOMODE);
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
    for (j=0;j<numstep;++j)
      rmsd_ave=(j*rmsd_ave+rmsd_trj[j])/(j+1);
    printf("cyc: %d  rmsd: %10.6lf\n",i+1,rmsd_ave);
    if (abs(rmsd_ave-rmsd_ave_old) < epsilon && i > 0)
      flago='o';
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
