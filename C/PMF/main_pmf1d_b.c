
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PMF.h"
#include "PT.h"
#include "EF.h"

double *io_scandcoldata(FILE *inputfile,int numi,int numcol,int col,int *numstep,double *data);

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,flag=0;

  int numatom,numstep;
  double *data;
  int numi=1,numcol=1,col=1;
  double max,min;
  int frame;
  double width=0.1;
  double *pmf;

  double T=0.0,beta=1.0;
  double k_B=1.98723e-3;
  double UNITT=418.4070;
  
  char *inputfilename,*inputfilename2,*outputfilename,*c;
  FILE *inputfile1,*inputfile2,*outputfile;

  char *line;
  size_t len=0;

  int d;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *progname,*logfilename;

  FILE *logfile;

  int opt_idx=1;

  struct option long_opt[] = {
    {"width",1,NULL,'w'},
    {"temp",1,NULL,'t'},
    {"beta",1,NULL,'b'},
    {"numi",1,NULL,'i'},
    {"numcol",1,NULL,'c'},
    {"col",1,NULL,'x'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((d=getopt_long(argc,argv,"hw:t:b:i:c:x:",long_opt,&opt_idx))!=-1) {
    switch(d) {
    case 'w':
      width=atof(optarg);
      break;
    case 't':
      T=atof(optarg);
      break;
    case 'b':
      beta=atof(optarg);
      break;
    case 'i':
      numi=atoi(optarg);
      break;
    case 'x':
      col=atoi(optarg);
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

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }  
  inputfilename  = *argv;
  outputfilename = *++argv;

  numstep=0;
  data=(double *)gcemalloc(sizeof(double)*2);
  inputfile1=efopen(inputfilename,"r");
  data=io_scandcoldata(inputfile1,numi,numcol,col,&numstep,data);
  fclose(inputfile1);

  pmf=pmf_1dmap(data,numstep,width,&max,&min,&frame);

  if (T>0) beta=1.0/(k_B*T); 
 
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<frame;++i) fprintf(outputfile,"%e %e \n",width*i+min,1.0/beta*pmf[i]);
  fclose(outputfile);
  
  return 0;
}

double *io_scandcoldata(FILE *inputfile,int numi,int numcol,int col,int *numstep,double *data){
  int i,j,k;
  double f;

  char *line;
  size_t len=0;

  for (i=0;i<numi;++i)
    getline(&line,&len,inputfile);

  for (i=(*numstep);;++i) {
    for (j=0;j<numcol;++j) {
      if (fscanf(inputfile,"%lf",&f)!=-1) {
	if (j==col-1) {
	  data=(double *)gcerealloc(data,sizeof(double)*(i+1));
	  data[i]=f;
	}
      }
      else {
	*numstep=i-1;
	return data;
      }
    }
  }
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilenam outputfilename\n",progname);
}
