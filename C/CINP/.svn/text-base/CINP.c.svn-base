
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "CINP.h"
#include "EF.h"

int compare_int(const void* a, const void* b);

int Create_inpcrd_from_trj(char *outputfilenamebase,FILE *inputfile,int numsnap, int numinpcrd, int numatom,int *num, int amberflag) {
  int i,j;
  double *crd;
  
  char outputfilename[100];
  FILE *outputfile;
  
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  j=0;
  for (i=0;i<numsnap;++i) {
    if (i!=num[j])
      io_dismisstrj(inputfile,1,numatom);
    else {
      io_scanconf(inputfile,numatom,crd,'x');
      sprintf(outputfilename,"%s_%d",outputfilenamebase,j+1);
      outputfile=efopen(outputfilename,"w");
      if (amberflag==ON)
	io_outputconf_Amberform(outputfile,numatom,crd);
      else
	io_outputconf(outputfile,numatom,crd);
      fclose(outputfile);
      ++j;
    }
  }
  
  
  return 0;
}

int Create_inpcrd_from_trj_Amber(char *outputfilenamebase,FILE *inputfile,int numsnap, int numinpcrd, int numatom,int *num, int amberflag) {
  int i,j;
  double *crd;
  
  char outputfilename[100];
  FILE *outputfile;
  char *line;
  size_t len=0;

  getline(&line,&len,inputfile);
  
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  j=0;
  for (i=0;i<numsnap;++i) {
    if (i!=num[j])
      io_dismisstrj(inputfile,1,numatom);
    else {
      io_scanconf(inputfile,numatom,crd,'x');
      sprintf(outputfilename,"%s_%d",outputfilenamebase,j+1);
      outputfile=efopen(outputfilename,"w");
      if (amberflag==ON)
	io_outputconf_Amberform(outputfile,numatom,crd);
      else
	io_outputconf(outputfile,numatom,crd);
      fclose(outputfile);
      ++j;
    }
  }
  
  
  return 0;
}


int Create_inpcrd_from_trj_with_ene_check(char *outputfilenamebase, char *inputfilename,int numsnap, int numinpcrd, int numatom,int *num, int amberflag, double *ene, double checkv) {
  int i,j,k,l;
  double *crd;
  double cvalue;
  double *qsene;
  
  char outputfilename[100];
  FILE *inputfile,*outputfile,*log1,*log2;

  qsene=(double *)gcemalloc(sizeof(double)*numsnap);
  memcpy(qsene,ene,sizeof(double)*numsnap);
  qsort((void *)qsene,numsnap,sizeof(double),compare_int);
  if (checkv < 0.0 || checkv > 1.0){
    printf("error:checkv must be [0.0:1.0]\n");
    exit(1);
  }
  cvalue=qsene[(int)(numsnap*checkv)-1];

  log1=efopen("log_ene.txt","w");
  log2=efopen("log_crd.txt","w");
  crd=(double *)gcemalloc(sizeof(double)*numatom*3);
  j=0;
  while (j<numinpcrd) {
    l=0;
    inputfile=efopen(inputfilename,"r");  
    for (i=0;i<numsnap;++i) {
      if (num[l]>=numsnap) exit(1);
      if (i!=num[l] || ene[i]>=cvalue )
	io_dismisstrj(inputfile,1,numatom);
      else {
	io_scanconf(inputfile,numatom,crd,'x');
	sprintf(outputfilename,"%s_%d",outputfilenamebase,j+1);
	outputfile=efopen(outputfilename,"w");
	if (amberflag==ON) {
	  io_outputconf_Amberform(outputfile,numatom,crd);
	  io_outputconf(log2,numatom,crd);
	}
	else {
	  io_outputconf(outputfile,numatom,crd);
	  io_outputconf(log2,numatom,crd);
	}
	fprintf(log1,"%d %e\n",i,ene[i]);
	fclose(outputfile);
	++j;++l;
      }
    }
    for (k=0;k<numinpcrd;++k) num[k]+=1;
    fclose(inputfile);
  }

  fclose(log1);
  fclose(log2);
  
  return 0;
}

int compare_int(const void* a, const void* b){
  if (*(double*)a < *(double*)b)
    return -1;
  else if (*(double*)a == *(double*)b)
    return 0;
  else 
    return 1;
}
