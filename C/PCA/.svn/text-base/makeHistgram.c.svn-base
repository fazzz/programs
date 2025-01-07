#include <stdio.h>
#include <stdlib.h>

#include "Histgram.h"

int makeHistgram(double *data, int datanum, double width, double max, double min, struct histgram hist){
  int i,j;
  int size;

  // finding minimam value
  min = data[0];
  max = data[0];
  for (i=0;i<datanum;++i){
    if (min>data[i]){
      min=data[i];
    }
    if (max<data[i]){
      max=data[i];
    }
  }
  size=(max-min)/width+1;
  hist->width=malloc(sizeof(double)*size);
  hist->height=malloc(sizeof(int)*size);

  for (j=0;j<size;++j) {
    hist.width[j]=min+width*j;
  }

  for (i=0;i<datanum;++i) {
    for (j=0;j<size;++j) {
      if (min+width*j<=data[i] && min+width*(j+1)>data[i]){
	hist.width[j]=min+width*j;
	hist.height[j]+=1;
      }
    }
  }

}
