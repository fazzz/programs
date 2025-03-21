
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "HIST.h"
#include "EF.h"
#include "PMF.h"

double kb=1.98723e-3*4.18407*100.0;

double *pmf_2dmap(double *data, int numdata,double width,double *maxx,double *maxy,double *minx,double *miny,int *framex,int *framey ){
  int i,j;
  double *pmf,pmf_min=0.0;

  pmf=hist_mk2dhist(data,width,width,numdata,maxx,maxy,minx,miny,framex,framey,1);

  for (i=0;i<=*framex;++i) {
    for (j=0;j<*framey;++j) {
      if ( pmf_min == 0.0 && pmf[i*(*framey)+j] > 0.0 )
	pmf_min=pmf[i*(*framey)+j];
      if ( pmf_min > pmf[i*(*framey)+j] && pmf[i*(*framey)+j] > 0.0 )
	pmf_min=pmf[i*(*framey)+j];
    }
  }

  for (i=0;i<=*framex;++i)
    for (j=0;j<*framey;++j)
      if (pmf[i*(*framey)+j]!=0.0)
	pmf[i*(*framey)+j]=-log(pmf[i*(*framey)+j])+log(pmf_min);

  return pmf;
} 

double *pmf_2dmap_rew(double *data,double *w, int numdata,double width,double *maxx,double *maxy,double *minx,double *miny,int *framex,int *framey,double temp){
  int i,j;
  double *pmf,pmf_min=0.0;
  double *hist,sum;

  *maxx=data[0];*minx=data[0];
  *maxy=data[1];*miny=data[1];
  for (i=0;i<numdata;++i) {
    if (*maxx < data[i*2]) *maxx=data[i*2];
    if (*minx > data[i*2]) *minx=data[i*2];
    if (*maxy < data[i*2+1]) *maxy=data[i*2+1];
    if (*miny > data[i*2+1]) *miny=data[i*2+1];
  }

  *framex=(int)((*maxx-*minx)/width)+1;
  *framey=(int)((*maxy-*miny)/width)+1;

  hist=(double *)gcemalloc(sizeof(double)*(*framex)*(*framey));
  pmf=(double *)gcemalloc(sizeof(double)*(*framex)*(*framey));
  for (i=0;i<*framex;++i)
    for (j=0;j<*framey;++j)
      hist[i*(*framey)+j]=0.0;

  for (i=0;i<numdata;++i)
    if (data[2*i] <= *maxx && data[2*i+1] <= *maxy)
      hist[((int)((data[i*2]-*minx)/width))*(*framey)+((int)((data[i*2+1]-*miny)/width))]+=exp(w[i]/kb/temp);

  sum=0.0;
  for (i=0;i<*framex;++i)
    for (j=0;j<*framey;++j)
      sum+=hist[i*(*framey)+j];
  for (i=0;i<*framex;++i)
    for (j=0;j<*framey;++j)
      pmf[i*(*framey)+j]=hist[i*(*framey)+j]/sum;

  for (i=0;i<=*framex;++i) {
    for (j=0;j<*framey;++j) {
      if ( pmf_min == 0.0 && pmf[i*(*framey)+j] > 0.0 )
	pmf_min=pmf[i*(*framey)+j];
      if ( pmf_min > pmf[i*(*framey)+j] && pmf[i*(*framey)+j] > 0.0 )
	pmf_min=pmf[i*(*framey)+j];
    }
  }

  for (i=0;i<=*framex;++i)
    for (j=0;j<*framey;++j)
      if (pmf[i*(*framey)+j]!=0.0)
	pmf[i*(*framey)+j]=-log(pmf[i*(*framey)+j])+log(pmf_min);

  return pmf;
} 

double *pmf_1dmap(double *data, int numdata,double width,double *max,double *min,int *frame){
  int i,j;
  double *pmf,pmf_min=0.0;

  pmf=hist_mkhist(data,width,numdata,max,min,frame,1);

  for (i=0;i<=*frame;++i) {
    if ( pmf_min == 0.0 && pmf[i] > 0.0 )  pmf_min=pmf[i];
    if ( pmf_min > pmf[i] && pmf[i] > 0.0 ) pmf_min=pmf[i];
  }

  for (i=0;i<=*frame;++i) if (pmf[i]!=0.0) pmf[i]=-log(pmf[i])+log(pmf_min);

  return pmf;
} 

double *pmf_2dmap_wmaxmin(double *data, int numdata,double widthx,double widthy,double maxx,double maxy,double minx,double miny,int *framex,int *framey ){
  int i,j;
  double *pmf,pmf_min=0.0;

  pmf=hist_mk2dhist_wmaxmin(data,widthx,widthy,numdata,maxx,maxy,minx,miny,framex,framey,1);

  /******************************************************************************/
  /* //  for (i=0;i<=framex;++i) {					        */
  /* for (i=0;i<=*framex;++i) {						        */
  /*   //    for (j=0;j<framey;++j) {					        */
  /*   for (j=0;j<*framey;++j) {					        */
  /*     if ( pmf_min == 0.0 && pmf[i*(*framey)+j] > 0.0 )		        */
  /* 	pmf_min=pmf[i*(*framey)+j];					        */
  /*     if ( pmf_min > pmf[i*(*framey)+j] && pmf[i*(*framey)+j] > 0.0 )        */
  /* 	pmf_min=pmf[i*(*framey)+j];					        */
  /*   }								        */
  /* }									        */
  /* 									        */
  /* //  for (i=0;i<=framex;++i)					        */
  /* for (i=0;i<=*framex;++i)						        */
  /*   //    for (j=0;j<framey;++j)					        */
  /*   for (j=0;j<*framey;++j)						        */
  /*     if (pmf[i*(*framey)+j]!=0.0)					        */
  /* 	pmf[i*(*framey)+j]=-log(pmf[i*(*framey)+j])+log(pmf_min);	        */
  /******************************************************************************/
  
  return pmf;
} 

double *pmf_2dmap_wmaxmin_rew(double *data, double *w, int numdata,double widthx,double widthy,double maxx,double maxy,double minx,double miny,int *framex,int *framey, double T ){
  int i,j;
  double *pmf,pmf_min=0.0;

  *framex=(int)((maxx-minx)/widthx)+1;
  *framey=(int)((maxy-miny)/widthy)+1;

  pmf=(double *)gcemalloc(sizeof(double)*(*framex)*(*framey));
  for (i=0;i<*framex;++i)
    for (j=0;j<*framey;++j)
      pmf[i*(*framey)+j]=0.0;

  for (i=0;i<numdata;++i)
    if (data[2*i] <= maxx && data[2*i+1] <= maxy)
      pmf[((int)((data[i*2]-minx)/widthx))*(*framey)+((int)((data[i*2+1]-miny)/widthy))]+=exp(w[i]/kb/T);

  for (i=0;i<*framex;++i)
    for (j=0;j<*framey;++j)
      pmf[i*(*framey)+j]/=numdata;
  
  return pmf;
} 
