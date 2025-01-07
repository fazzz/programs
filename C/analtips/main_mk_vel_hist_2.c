#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "netcdf_mineL.h"

#include "HIST.h"

#include "PTL.h"
#include "EF.h"
#include "IO.h"

#define ON 1
#define OFF 0

int usage(void);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int numatom,numspatom,numstep,initialstep=0;
  
  int normflag='o';
  int frame;
  double width=0.1;
  double maxx,minx,maxy,miny,maxz,minz;

  double *vel;
  double *datax,*datay,*dataz,*histx,*histy,*histz;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*outputfilename,*parmfilename,*histfilename;
  char *progname;
  FILE *inputfile,*parmfile,*outputfile,*histfile;

  while((c=getopt(argc,argv,"hn:w:i:"))!=-1) {
    switch(c) {
    case 'i':
      initialstep=atoi(optarg[0]);
    case 'n':
      normflag=optarg[0];
      break;
    case 'w':
      width=atof(optarg);
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

  if (argc < 6) {
    USAGE(progname);
    exit(1);
  }
  numspatom      = atoi(*argv);
  numstep        = atoi(*++argv);
  inputfilename  = *++argv;
  parmfilename   = *++argv;
  outputfilename = *++argv;
  histfilename   = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile);
  numatom=AP.NATOM;

  vel=(double *)gcemalloc(sizeof(double)*numatom*3);

  datax=(double *)gcemalloc(sizeof(double)*(numstep-initialstep));
  datay=(double *)gcemalloc(sizeof(double)*(numstep-initialstep));
  dataz=(double *)gcemalloc(sizeof(double)*(numstep-initialstep));

  inputfile=efopen(inputfilename,"r");
  outputfile=efopen(outputfilename,"w");
  //  getline(&line,&len,inputfile);

  for (i=0;i<numstep;++i) {
    for (j=0;j<numatom;++j) {
      for (k=0;k<3;++k) {
	fscanf(inputfile,"%lf ",&vel[j*3+k]);	
      }
      if ( i >= initialstep && j==numspatom-1) {
	datax[i-initialstep]=vel[j*3];
	datay[i-initialstep]=vel[j*3+1];
	dataz[i-initialstep]=vel[j*3+2];
	for (k=0;k<3;++k) {
	  fprintf(outputfile,"%lf ",vel[j*3+k]);
	}
	fprintf(outputfile,"\n");
      }
    }
  }
  fclose(inputfile);
  fclose(outputfile);

  histx=hist_mkhist(datax,width,(numstep-initialstep),&maxx,&minx,&frame,normflag);
  histy=hist_mkhist(datay,width,(numstep-initialstep),&maxy,&miny,&frame,normflag);
  histz=hist_mkhist(dataz,width,(numstep-initialstep),&maxz,&minz,&frame,normflag);

  histfile=efopen(histfilename,"w");
  for (i=0;i<=frame;++i) {
    fprintf(histfile,"%e %e  ",width*(i+0.5)+minx,histx[i]);
    fprintf(histfile,"%e %e  ",width*(i+0.5)+miny,histy[i]);
    fprintf(histfile,"%e %e\n",width*(i+0.5)+minz,histz[i]);
  }
  fclose(histfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h help]:\n");
  printf("%s inputfilename parmfilename outputfilename histfilename\n",progname);
}

 
