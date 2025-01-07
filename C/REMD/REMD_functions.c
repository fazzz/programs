
#define _GNU_SOURCE  

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "REMD_functions.h"

#define ON  0
#define OFF 1

int REMD_purmutation_func(int i){
  int m;

  m=index_parameters[i];

  return m; 
}

int REMD_purmutation_inverse_func(int m){
  int i;

  i=index_replicas[m];

  return i;     
}

int REMD_exchange_purmutation_funcs(int i,int j,int m, int n){
  int temp;

  if (i == -1 || j == -1 || m == -1 || n == -1) {

  }
  else {
    //  printf("i=%d j=%d m=%d n=%d index_parameters[i]=%d index_parameters[j]=%d\n",
    //	 i,j,m,n,index_parameters[i],index_parameters[j]);
    if (index_parameters[i]==m) index_parameters[i]=n;
    if (index_parameters[j]==n) index_parameters[j]=m;
    //  printf("index_parameters[i]=%d index_parameters[j]=%d\n",
    //	 index_parameters[i],index_parameters[j]);

    //  printf("index_replicas[m]=%d index_replicas[n]=%d\n",
    //	 index_replicas[m],index_replicas[n]);
    if (index_replicas[m]==i)  index_replicas[m]=j;
    if (index_replicas[n]==j)  index_replicas[n]=i;
    //  printf("index_replicas[m]=%d index_replicas[n]=%d\n",
    //	 index_replicas[m],index_replicas[n]);
  }
}

int REMD_ini_purmutation_funcs(int num){
  int i;

  for (i=0;i<num;++i) {
    index_replicas[i]=i;
    index_parameters[i]=i;
  }
}
