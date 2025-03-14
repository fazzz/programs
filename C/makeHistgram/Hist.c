#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "HIST.h"
#include "EF.h"

double *hist_mkhist(double *data, double width, int numdata, double *max, double *min, int *frame, int normflag) {
  int i;
  double *hist;

  *max=data[0];*min=data[0];
  for (i=0;i<numdata;++i) {
    if (*max < data[i]) *max=data[i];
    if (*min > data[i]) *min=data[i];
  }
  //  *frame=(int)((*max-*min)/width)+1;
  *frame=(int)round((*max-*min)/width)+1; // 2011-11-11
  hist=(double *)gcemalloc(sizeof(double)*(*frame));
  for (i=0;i<*frame;++i)
    hist[i]=0.0;

  for (i=0;i<numdata;++i)
    hist[((int)round((data[i]-*min)/width))]+=1.0; // 2011-11-11
  //  hist[((int)((data[i]-*min)/width))]+=1.0; 

  if (normflag==ON)
    for (i=0;i<*frame;++i)
      hist[i]/=numdata;

  return hist;

}


double *hist_mk2dhist(double *data, double widthx, double widthy , int numdata, double *maxx, double *maxy, double *minx, double *miny,int *framex, int *framey, int normflag) {
  int i,j;
  double *hist;

  *maxx=data[0];*minx=data[0];
  *maxy=data[1];*miny=data[1];
  for (i=0;i<numdata;++i) {
    if (*maxx < data[i*2]) *maxx=data[i*2];
    if (*minx > data[i*2]) *minx=data[i*2];
    if (*maxy < data[i*2+1]) *maxy=data[i*2+1];
    if (*miny > data[i*2+1]) *miny=data[i*2+1];
  }

  *framex=(int)((*maxx-*minx)/widthx)+1;
  *framey=(int)((*maxy-*miny)/widthy)+1;

  hist=(double *)gcemalloc(sizeof(double)*(*framex)*(*framey));
  for (i=0;i<*framex;++i)
    for (j=0;j<*framey;++j)
      hist[i*(*framey)+j]=0.0;

  for (i=0;i<numdata;++i)
    if (data[2*i] <= *maxx && data[2*i+1] <= *maxy)
      hist[((int)((data[i*2]-*minx)/widthx))*(*framey)+((int)((data[i*2+1]-*miny)/widthy))]+=1.0;

  if (normflag==ON)
    for (i=0;i<*framex;++i)
      for (j=0;j<*framey;++j)
	hist[i*(*framey)+j]/=numdata;

  return hist;
}

double *hist_mk2dhist_ts(double *data, double widthx, double widthy , int numdata, double *maxx, double *maxy, double *minx, double *miny,int *framex, int *framey) {
  int i,j,k;
  double *hist_ts;

  *maxx=data[0];*minx=data[0];
  *maxy=data[1];*miny=data[1];
  for (i=0;i<numdata;++i) {
    if (*maxx < data[i*2]) *maxx=data[i*2];
    if (*minx > data[i*2]) *minx=data[i*2];
    if (*maxy < data[i*2+1]) *maxy=data[i*2+1];
    if (*miny > data[i*2+1]) *miny=data[i*2+1];
  }

  *framex=(int)((*maxx-*minx)/widthx)+1;
  *framey=(int)((*maxy-*miny)/widthy)+1;

  hist_ts=(double *)gcemalloc(sizeof(double)*numdata*2);

  for (i=0;i<numdata;++i) {
    hist_ts[i*2]=((int)((data[i*2]-*minx)/widthx));
    hist_ts[i*2+1]=((int)((data[i*2+1]-*miny)/widthy));
  }

  return hist_ts;
}

double *hist_mkhist_wp(double *data, double width, int numdata, double *max, double *min, int *frame, int normflag, double *prob) {
  int i;
  double *hist;

  *max=data[0];*min=data[0];
  for (i=0;i<numdata;++i) {
    if (*max < data[i]) *max=data[i];
    if (*min > data[i]) *min=data[i];
  }
  *frame=(int)((*max-*min)/width)+1;
  hist=(double *)gcemalloc(sizeof(double)*(*frame));
  for (i=0;i<*frame;++i)
    hist[i]=0.0;

  for (i=0;i<numdata;++i)
    hist[((int)((data[i]-*min)/width))]+=prob[i];

  if (normflag==ON)
    for (i=0;i<*frame;++i)
      hist[i]/=numdata;

  return hist;

}

double *hist_mk2dhist_wmaxmin(double *data, double widthx, double widthy , int numdata, double maxx, double maxy, double minx, double miny,int *framex, int *framey, int normflag) {
  int i,j;
  double *hist;

  *framex=(int)((maxx-minx)/widthx)+1;
  *framey=(int)((maxy-miny)/widthy)+1;

  hist=(double *)gcemalloc(sizeof(double)*(*framex)*(*framey));
  for (i=0;i<*framex;++i)
    for (j=0;j<*framey;++j)
      hist[i*(*framey)+j]=0.0;

  for (i=0;i<numdata;++i)
    if (data[2*i] <= maxx && data[2*i+1] <= maxy)
      hist[((int)((data[i*2]-minx)/widthx))*(*framey)+((int)((data[i*2+1]-miny)/widthy))]+=1.0;

  if (normflag==ON)
    for (i=0;i<*framex;++i)
      for (j=0;j<*framey;++j)
	hist[i*(*framey)+j]/=numdata;

  return hist;
}
