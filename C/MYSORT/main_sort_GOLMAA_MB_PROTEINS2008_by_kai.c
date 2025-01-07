#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "PTL.h"
#include "FFL.h"
#include "EF.h"

#define ARRAY_SIZE(array) (sizeof(array)/sizeof(array[0]))

struct ene_kai {
  double ene;
  double kai;
};

typedef struct ene_kai EK;

int compare_ene_kai( const void *s1 , const void *s2 );

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,n;
  double f1,f2,da=0.1,kai,dkai=0.1;
  int minindex;
  double minene;
  struct ene_kai* ene_kai_data;


  double pi;

  char *inputfilename;
  char *outputfilename;

  FILE *inputfile;
  FILE *outputfile;

  FILE *logfile;

  char *progname;

  int opt_idx=1;
  int c;

  struct option long_opt[] = {
    {"da",1,NULL,'@'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"h@:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case '@':
      da=atof(optarg);
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  pi=acos(-1.0);

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 2) {
    USAGE(progname);
    exit(1);
  }
  inputfilename  = *argv;
  outputfilename = *++argv;

  ene_kai_data=(struct ene_kai *)gcemalloc(sizeof(struct ene_kai)*1);

  inputfile=efopen(inputfilename,"r");
  for (i=0;(c=fscanf(inputfile,"%lf",&f1))!=-1;++i) {
    fscanf(inputfile,"%lf",&f2);
    ene_kai_data=(struct ene_kai *)gcerealloc(ene_kai_data,sizeof(struct ene_kai)*(i+1));
    ene_kai_data[i].ene=f1;
    ene_kai_data[i].kai=f2;
  }
  fclose(inputfile);
  n=i;

  qsort(ene_kai_data,n,sizeof(struct ene_kai),compare_ene_kai);

  outputfile=efopen(outputfilename,"w");
  kai=ene_kai_data[0].kai;
  minindex=0;
  for (i=0;i<n;) {
    minene=ene_kai_data[i].ene;
    for (;ene_kai_data[i].kai<kai;++i) {
      if (minene > ene_kai_data[i].ene) {
	minindex=i;
      }
    }
    fprintf(outputfile,"%e %e \n",ene_kai_data[minindex].kai,ene_kai_data[minindex].ene);
    kai+=dkai;
    ++i;
  }
  fclose(outputfile);

}

void USAGE(char *progname) {
   printf("-h -- help\n");
   printf("USAGE: %s inputfilename outputfilename\n", progname);
 }

int compare_ene_kai( const void *s1 , const void *s2 ) { 
  EK *a1,*a2;

  a1=(EK *)s1;
  a2=(EK *)s2;

  if( a1->kai < a2->kai ){ 
    return -1; 
  } else if( a1->kai == a2->kai ){ 
    return 0; 
  } else { 
    return 1; 
  } 
}

