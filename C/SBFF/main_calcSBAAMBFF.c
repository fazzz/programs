
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "EF.h"
#include "IO.h"
#include "SBFF.h"
#include "PT.h"
#include "netcdf_mine.h"

#define ON 1
#define OFF 0

#define kb 1.98723e-3

void USAGE(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int d;
  double f;
  int numarg;

  int numatom;

  double pi;

  double delta=1.0,deltaV;
  double p_t;

  double *crd1,*crd2;
  struct potential_SBAA ene1,ene2;
  struct potential_SBAAMB ene;

  char *progname;
  char *inifile1name,*inifile2name,*parmtopname,*outfilename;
  FILE *inifile1,*inifile2,*parmtop,*outfile;

  struct my_netcdf_out_id_MCD nc_id_MCD;

  char *line,*dummy;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  progname=argv[0];
  while((c=getopt(argc,argv,"h"))!=-1) {
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

  numarg=4;

  if (argc < numarg) {
    USAGE(progname);
    exit(1);
  }
  inifile1name  = *argv;
  inifile2name  = *++argv;
  parmtopname  = *++argv;
  outfilename  = *++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtop(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  crd1=(double *)gcemalloc(sizeof(double )*numatom*3);
  crd2=(double *)gcemalloc(sizeof(double )*numatom*3);

  inifile1=efopen(inifile1name,"r");
  io_scanconf_Amber_ini(inifile1,numatom,crd1);
  fclose(inifile1);

  inifile2=efopen(inifile2name,"r");
  io_scanconf_Amber_ini(inifile2,numatom,crd2);
  fclose(inifile2);

  SBAAff_set_parameters_default(&ene1);
  SBAAff_set_parameters_default(&ene2);

  SBAAMBff_set_calcff(&ene1,&ene2,&deltaV,crd1,crd2,numatom);

  outfile=efopen(outfilename,"w");

  p_t=SBAAMBff_calcff(crd1,numatom,&ene1,&ene2,delta,deltaV);
  fprintf(outfile,"%e\n",p_t);

  p_t=SBAAMBff_calcff(crd2,numatom,&ene1,&ene2,delta,deltaV);
  fprintf(outfile,"%e\n",p_t);

  fclose(outfile);

}
 
void USAGE(char *progname) {
  printf("[-h] -- help\n");
  printf("USAGE: %s inifile1name inifile2name parmtopname outfilename\n", progname);
}

