
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "mymath.h"
#include "PT.h"
#include "IO.h"
#include "EF.h"

int usage(void);

int main(int argc, char *argv[]) {
  int i,j;
  int flag,flag2,flag3,c;
  int nums,numv; 
  double *timeseries,*timeseriesnorm;

  char *line;
  size_t len=0;

  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *inputfilename1,*inputfilename2,*outputfilename;
  FILE *inputfile1,*inputfile2,*outputfile;

  while((c=getopt(argc,argv,"LFM"))!=-1) {
    switch(c) {
    case 'L':
      flag='L';
      break;
    case 'F':
      flag2='F';
      break;
    case 'M':
      flag3='M';
      break;
    default:
      usage();
      exit(1);
    }
  }

  argc -= optind;
  argv += optind;

  if (argc < 4 ) {
    printf("USAGE: %s nums nuumv inputfilename1(data) outputfilename(data_norm)\n",argv[0]);
    exit(1);
  }
  nums=atoi(*argv);
  numv=atoi(*++argv);
  inputfilename1 = *++argv;
  outputfilename = *++argv;
  
  inputfile1=efopen(inputfilename1,"r");
  if (flag=='L')
    getline(&line,&len,inputfile1);
  timeseries=(double *)gcemalloc(sizeof(double)*numv*nums);
  timeseriesnorm=(double *)gcemalloc(sizeof(double)*numv*nums);
  if (flag3=='M')
    io_scantimeseries(inputfile1,nums,numv,timeseries,'c');
  else
    io_scantimeseries(inputfile1,nums,numv,timeseries,'x');
  fclose(inputfile1);

  ts_normalize(timeseries,nums,numv,timeseriesnorm);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<nums;++i) {
    if (flag2=='F')
      fprintf(outputfile,"%d ",i+1);
    for (j=0;j<numv;++j)
      fprintf(outputfile,"%12.8e ",timeseriesnorm[i*numv+j]);
    fprintf(outputfile,"\n ");
  }
  fclose(outputfile);

  return 0;
}

int usage(void) {
  printf("USAGE:norm \n");
  printf("-L \n");
  printf("nums nuumv inputfilename1(data) outputfilename(data_norm)\n");
}
