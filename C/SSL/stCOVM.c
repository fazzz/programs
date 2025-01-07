
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "efunc.h"
#include "SSL.h"

double ssl_covm(double *timeseries, int nums, int numv, double *COVM){
  int i,j,k;

  COVM=ecalloc(numv*numv,sizeof(double));
  for (i=0;i<nums;++i)
    for (j=0;j<numv;++j)
      for (k=j;k<numv;++k)
	COVM[i*numv+j]=(i*COVM[i*numv+j]+timeseries[i*numv*numv+j*numv+k]*timeseries[i*numv*numv+j*numv+k])/(i+1);

  return 0.0;
}

double ssl_ave(double *timeseries, int nums, int numv, double *AVE){
  int i,j;

  AVE=ecalloc(numv,sizeof(double));
  for (i=0;i<nums;++i)
    for (j=0;j<numv;++j)
      AVE[i]=i*AVE[i]+timeseries[i*numv+j]/(i+1);

  return 0.0;
}


double normalize(double *timeseries, int nums, int numv) {

  ssl_ave(timeseries,nums,numv,ave);
  for (i=0;i<numv;++i)
    for (i=0;i<numv;++i)
      timeseries_n[i]=timeseries[i*nums+j]-ave[i];


}

