
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "SEQUENCE.h"
#include "EF.h"

double *geome_seq_wmin_max(double min, double max, double common ){
  int i=0;
  double *seq;

  seq=(double *)gcemalloc(sizeof(double)*1);

  while (seq[i]<max) {
    ++i;
    seq=(double *)gcerealloc(seq,sizeof(double)*(i+1));
    seq[i]=seq[i-1]*common;
  }

  return seq;
}

double *geome_seq_wmin_num(double min, int num, double common ){
  int i=0;
  double *seq;

  seq=(double *)gcemalloc(sizeof(double)*num);

  seq[i]=min;
  for (i=1;i<num;++i) {
    seq[i]=seq[i-1]*common;
  }

  return seq;
}

double *arith_seq_wmin_max(double min, double max, double common ){
  int i=0;
  double *seq;

  seq=(double *)gcemalloc(sizeof(double)*1);

  while (seq[i]<max) {
    ++i;
    seq=(double *)gcerealloc(seq,sizeof(double)*(i+1));
    seq[i]=seq[i-1]+common;
  }

  return seq;
}

double *arith_seq_wmin_num(double min, int num, double common ){
  int i=0;
  double *seq;

  seq=(double *)gcemalloc(sizeof(double)*num);

  seq[i]=min;
  for (i=1;i<num;++i) {
    seq[i]=seq[i-1]+common;
  }

  return seq;
}
