
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "EF.h"
#include "PROTOPO.h"

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,n;
  int aromaflag=ON;

  char *name_atom_list;

  int *bp;
  int **bp_f;
 
  int **bpairs,*numb;
  int **matrix,numatom;
  
  int **pairex1_5,**pair1_4;
  int *num1_5, *num14;

  char *line;
  size_t len=0;
  
  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;
  
  char *progname;
  char *inputfilename;
  char *outputfilename;
  
  FILE *inputfile;
  FILE *outputfile;
  
  while((c=getopt(argc,argv,"ha"))!=-1) {
    switch(c) {
    case 'a':
      aromaflag=OFF;
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
  inputfilename = *argv;
  outputfilename = *++argv;

  inputfile=efopen(inputfilename,"r");
  fscanf(inputfile,"%d",&numatom);
  name_atom_list=(char *)gcemalloc(sizeof(char)*numatom*4);
  getline(&line,&len,inputfile);
  for (i=0;i<numatom;++i) {
    getline(&line,&len,inputfile);
    for (j=0;j<4;++j) name_atom_list[i*4+j]=line[j];
  }
  fclose(inputfile);

  bp=(int *)gcemalloc(sizeof(int)*(numatom-1)*2);
  bp_f=(int **)gcemalloc(sizeof(int *)*numatom);

  numb=(int *)gcemalloc(sizeof(int)*numatom);
  make_bp(numatom,bp,bp_f,numb,name_atom_list,aromaflag);
  
  for (i=0;i<numatom;++i) {
    for (j=0;j<4;++j) printf("%c",name_atom_list[i*4+j]);
    printf("(%3d)-",i);
    for (j=0;j<numb[i];++j) {
      n=bp_f[i][j];
      for (k=0;k<4;++k) printf("%c",name_atom_list[n*4+k]);
      printf("(%3d) ",n);
    }
    printf("\n");
  }

  matrix=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) matrix[i]=(int *)gcemalloc(sizeof(int)*numatom);

  make_nb_matrix(bp_f,numb,3,matrix,numatom);

  pair1_4=(int **)gcemalloc(sizeof(int *)*numatom);
  pairex1_5=(int **)gcemalloc(sizeof(int *)*numatom);
  num14=(int *)gcemalloc(sizeof(int)*numatom);
  num1_5=(int *)gcemalloc(sizeof(int)*numatom);

  set_nb_pairs(matrix,numatom,pairex1_5,pair1_4,num1_5,num14);

  for (i=0;i<numatom;++i) {
    //    for (j=0;j<4;++j) printf("%c",name_atom_list[i*4+j]);
    //    printf("--");
    //    for (j=0;j<num14[i];++j) 
    //      for (k=0;k<4;++k) printf("%c",name_atom_list[(pair1_4[i][j])*4+k]);
    //    printf("\n");
  }

  for (i=0;i<numatom;++i) {
    //    for (j=0;j<4;++j) printf("%c",name_atom_list[i*4+j]);
    //    printf("--");
    //    for (j=0;j<num1_5[i];++j) 
    //      for (k=0;k<4;++k) printf("%c",name_atom_list[(pairex1_5[i][j])*4+k]);
    //    printf("\n");
  }

  for (i=0;i<numatom;++i) {
    //    for (j=0;j<4;++j) printf("%c",name_atom_list[i*4+j]);
    for (j=0;j<num14[i];++j)
      printf("MA=%2d MB=%2d\n",i+1,pair1_4[i][j]+1);
  }
  
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numatom;++i) {
    for (j=0;j<numatom;++j) 
      fprintf(outputfile,"%4d",matrix[i][j]);
    fprintf(outputfile,"\n");
  }
  fclose(outputfile);
  
  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s inputfilename outputfilename \n",progname);
}
