
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

  double *crd;
  struct potential_SBAA ene,ene_trial;

  char *progname;
  char *inifilename,*parmtopname,*outfilename;
  FILE *inifile,*parmtop,*outfile;

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

  numarg=3;

  if (argc < numarg) {
    USAGE(progname);
    exit(1);
  }
  inifilename  = *argv;
  parmtopname  = *++argv;
  outfilename  = *++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtop(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;

  crd=(double *)gcemalloc(sizeof(double )*numatom*3);

  inifile=efopen(inifilename,"r");
  io_scanconf_Amber_ini(inifile,numatom,crd);
  fclose(inifile);

  SBAAff_set_parameters_default(&ene);

  SBAAff_set_calcff(&ene,crd,numatom);

  SBAAff_check_parameters(&ene,numatom);

  SBAAff_calcff(crd,numatom,&ene);

  outfile=efopen(outfilename,"w");
  fprintf(outfile,"%e\n",ene.p_t);
  fprintf(outfile,"%e\n",ene.p_cnb_t);
  fprintf(outfile,"%e\n",ene.p_nnb_t);
  fprintf(outfile,"%e\n",ene.p_d_t);
  fprintf(outfile,"%e\n",ene.p_a_t);
  fprintf(outfile,"%e\n",ene.p_b_t);
  for (i=0;i<numatom;++i) {
    fprintf(outfile,"%e\n",ene.p_cnb[i]);
  }
  for (i=0;i<numatom;++i) {
    fprintf(outfile,"%e\n",ene.p_nnb[i]);
  }
  for (i=0;i<AP.MBONA;++i) {
    fprintf(outfile,"%e\n",ene.p_b[i]);
  }
  for (i=0;i<AP.MTHETA;++i) {
    fprintf(outfile,"%e\n",ene.p_a[i]);
  }
  for (i=0;i<AP.MPHIA;++i) {
    fprintf(outfile,"%e\n",ene.p_d[i]);
  }
  fclose(outfile);

}
 
void USAGE(char *progname) {
  printf("[-h] -- help\n");
  printf("USAGE: %s inifilename parmtopname outfilename\n", progname);
}

