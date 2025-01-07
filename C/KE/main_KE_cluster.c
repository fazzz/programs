
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "KE.h"
#include "PT.h"
#include "EF.h"
#include "IO.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m,d1,d2;
  int numstep,numclut;
  double ke;

  int *nainiclt,**nafinclt;

  double dt=0.001;
  struct clust clt;
  double *crd1,*crd2;
  double *crdclt1,*crdclt2;
  double *mass,*massclt;
  int numatom;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *progname;
  char *inputfilename,*parmtopname,*clustfilename,*massfilename;
  char *outputfilename;

  FILE *inputfile,*parmtop,*clustfile,*massfile;
  FILE *outputfile;

  while((c=getopt(argc,argv,"hd:"))!=-1) {
    switch(c) {
    case 'd':
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

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename = *argv;
  clustfilename = *++argv;
  massfilename  = *++argv;
  outputfilename = *++argv;

  /************************************/
  /* parmtop=efopen(parmtopname,"r"); */
  /* readParmtop(parmtop);	      */
  /* fclose(parmtop);		      */
  /************************************/

  clustfile=efopen(clustfilename,"r");
  io_scanclust(clustfile,&clt);
  fclose(clustfile);
  numclut=clt.numclt;
  nainiclt=(int)gcemalloc(sizeof(int)*numclut);
  nafinclt=(int *)gcemalloc(sizeof(int *)*numclut);
  for (i=0;i<numclut;++i) nafinclt[i]=(int)gcemalloc(sizeof(int)*clt.nbranch[i]);
  numatom=0;
  for (i=0;i<numclut;++i) numatom+=clt.natom[i];

  massfile=efopen(massfilename,"r");
  mass=(double *)gcemalloc(sizeof(double)*numatom);
  massclt=(double *)gcemalloc(sizeof(double)*numatom);
  for (i=0;i<numatom;++i) fscanf(massfile,"%lf",&mass[i]);
  fclose(massfile);

  crd1=(double *)gcemalloc(sizeof(double)*numatom*3);
  crd2=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdclt1=(double *)gcemalloc(sizeof(double)*numatom*3);
  crdclt2=(double *)gcemalloc(sizeof(double)*numatom*3);


  for (i=0;i<numclut;++i) nainiclt[i]=clt.nainiclt[i];
  for (i=0;i<numclut;++i) {
    for (j=0;j<clt.nbranch[i];++j)
      nafinclt[i][j]=clt.nafinclt[i][j];
  }

  inputfile=efopen(inputfilename,"r");
  outputfile=efopen(outputfilename,"w");
  for (i=0;d1!=-1 && d2!=-1;++i) {
    d1=io_scanconfwj(inputfile,numatom,crd1,'x');
    d2=io_scanconfwj(inputfile,numatom,crd2,'x');

    for (j=0;j<numclut;++j) {
      m=0;
      for (k=m;k<nafinclt[j][0];++k) {
	massclt[m]=mass[k];
	for (l=0;l<3;++l) {
	  crdclt1[m*3+l]=crd1[k*3+l];
	  crdclt2[m*3+l]=crd2[k*3+l];
	}
	++m;
      }
      ke=kieneftrj(crdclt1,crdclt2,massclt,dt,m);
      fprintf(outputfile,"%e ",ke/numatom/3/kbkcl);
    }
    fprintf(outputfile,"\n");
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s inputfilename clustfilename outputfilename \n",progname);
}
