
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "NC.h"
#include "PTL.h"
#include "EF.h"

#include "netcdf_mine.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d,na,ii,jj;
  int numstep,HMODE=EXC;
  double QCA;

  double *crd,*crdref;
  int numatom,numres;

  double criteria=criteria_NC;
  int numncaa,numncres;
  int  *index_natatt,**ncmap,**ncmapres;
  double atom1[3],atom2[3];

  double *cradii_ca;
  double vec[3];
  double len;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *progname;
  char *crdfilename,*parmtopname,*crdreffilename;
  char *outputfilename;

  FILE *crdfile,*parmtop,*crdreffile;
  FILE *outputfile;

  while((c=getopt(argc,argv,"hc:"))!=-1) {
    switch(c) {
    case 'c':
      criteria=atof(argv);
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
  crdfilename = *argv;
  crdreffilename = *++argv;
  parmtopname = *++argv;
  outputfilename = *++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;
  numres=AP.NRES;

  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdref=(double *)gcemalloc(sizeof(double)*numatom*3);

  crdreffile=efopen(crdreffilename,"r");
  io_scanconf_Amber_ini(crdreffile,numatom,crdref);
  fclose(crdreffile);

  crdfile=efopen(crdfilename,"r");
  io_scanconf_Amber_ini(crdfile,numatom,crd);
  fclose(crdfile);

  ncmap=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) ncmap[i]=(int *)gcemalloc(sizeof(int)*numatom);
  ncmapres=(int **)gcemalloc(sizeof(int *)*numres);
  for (i=0;i<numres;++i) ncmapres[i]=(int *)gcemalloc(sizeof(int)*numres);

  index_natatt=make_native_contact_list_aa_3(&numncaa,&numncres,crdref,numatom,numres,criteria,ncmap,ncmapres,EXC);
  cradii_ca=(double *)gcemalloc(sizeof(double)*numncaa);

  na=0;
  for (i=0;i<numres;++i) {
    for (j=i+1;j<numres;++j) {
      if (ncmapres[i][j]==0 ) {
	ii=PTL_canum_fromresnum(i);
	jj=PTL_canum_fromresnum(j);
	len = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crdref[ii*3+k]-crdref[jj*3+k];
	  len += vec[k]*vec[k];
	}
	len = sqrt(len);
	cradii_ca[na]=len;
	printf("%d-%d %d-%d %lf\n",i,j,ii,jj,len);
	for (k=0;k<2;++k) printf("%c",AP.IGRAPH[ii][k]);
	printf(" ");
	for (k=0;k<2;++k) printf("%c",AP.IGRAPH[jj][k]);
	printf("\n ");
	++na;
      }
    }
  }

  outputfile=efopen(outputfilename,"w");
  QCA=count_native_contact_ca(numncres,crd,numres,ncmapres,cradii_ca);
  fprintf(outputfile,"%e\n",QCA);
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename reffilename parmfilename outputfilename \n",progname);
}

