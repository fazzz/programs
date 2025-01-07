
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "netcdf_mine.h"

#include "DHAD.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

int usage(void);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int numatom,numstep,*numdihed,numd;
  int modeflag='P';
  int **adpairs;

  double *crd,*crdref;
  double dhadis;
  int *list;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;

  char *inputfilename,*outputfilename,*parmfilename,*crdreffilename;
  char *progname;
  FILE *inputfile,*outputfile,*parmfile,*crdreffile;

  while((c=getopt(argc,argv,"hPOK"))!=-1) {
    switch(c) {
    case 'P':
      modeflag='P';
      break;
    case 'O':
      modeflag='O';
      break;
    case 'K':
      modeflag='K';
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  parmfilename   = *++argv;
  crdreffilename = *++argv;
  outputfilename = *++argv;

  numstep=mync_get_present_step_MCD(inputfilename,&nc_id_MCD);

  parmfile=efopen(parmfilename,"r");
  readParmtop(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdref=(double *)gcemalloc(sizeof(double)*numatom*3);

  crdreffile=efopen(crdreffilename,"r");
  io_scanconf_Amber_ini(crdreffile,numatom,crdref);
  //  io_scanconf_Amber_rst(crdreffile,crdref);
  fclose(crdreffile);

  adpairs=(int **)gcemalloc(sizeof(int *)*5);
  adpairs[0]=(int *)gcemalloc(sizeof(int)*4); // PHI,PSI
  adpairs[1]=(int *)gcemalloc(sizeof(int)*4); // OMEGA
  adpairs[2]=(int *)gcemalloc(sizeof(int)*4); // KI
  adpairs[3]=(int *)gcemalloc(sizeof(int)*8); // ACE
  adpairs[4]=(int *)gcemalloc(sizeof(int)*8); // NME
  numdihed=(int *)gcemalloc(sizeof(int)*5);  
  readdihedpairs(adpairs,numdihed);

  if (modeflag=='P') {
    list=(int *)gcemalloc(sizeof(int)*numdihed[0]*4);
    for (i=0;i<numdihed[0]*4;++i)   list[i]=adpairs[0][i];
    numd=numdihed[0];
  }
  else if (modeflag=='O') {
    list=(int *)gcemalloc(sizeof(int)*(numdihed[0]+numdihed[1])*4);
    for (i=0;i<numdihed[0]*4;++i)  list[i]=adpairs[0][i];
    for (i=0;i<numdihed[1]*4;++i)  list[numdihed[0]+i]=adpairs[1][i];
    numd=numdihed[0]+numdihed[1];
  }
  else {
    list=(int *)gcemalloc(sizeof(int)*(numdihed[0]+numdihed[1]+numdihed[2]+numdihed[3])*4);
    for (i=0;i<numdihed[0]*4;++i)  list[i]=adpairs[0][i];
    for (i=0;i<numdihed[1]*4;++i)  list[numdihed[0]+i]=adpairs[1][i];
    for (i=0;i<numdihed[2]*4;++i)  list[numdihed[1]+i]=adpairs[2][i];
    for (i=0;i<numdihed[3]*4;++i)  list[numdihed[2]+i]=adpairs[3][i];
    numd=numdihed[0]+numdihed[1]+numdihed[2]+numdihed[3];
  }
  
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numstep;++i) {
    mync_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);
    for (j=0;j<numatom;++j)
      for (k=0;k<3;++k) 
	crd[j*3+k]=crd_nc[j][k];
    dhadis=dhad(numd,list,crd,crdref);
    fprintf(outputfile,"%e\n",dhadis);
  }

  fclose(outputfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s inputfilename parmfilename crdreffilename outputfilename \n",progname);
}
