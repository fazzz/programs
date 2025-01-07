
#include <math.h>
#include "mymath.h"

double Fermi(double x) {
  return 1.0/(1.0+exp(x));
}

double calc_ave(int num, double *data) {
  int i;
  double ave=0.0;
  for (i=0;i<num;++i)
    ave=(i*ave+data[i])/(i+1);
  return ave;
}

double calc_msq(int num, double *data) {
  int i;
  double msq=0.0;
  for (i=0;i<num;++i)
    msq=(i*msq+data[i]*data[i])/(i+1);
  return msq;
}

double calc_var(int num, double *data) {
  return calc_msq(num,data)-calc_ave(num,data)*calc_ave(num,data);
}

double MYMATH_ave(double *data,int numstep) {
  int i;
  double ave=0.0;
  for (i=0;i<numstep;++i) ave+=data[i];
  ave=ave/numstep;

  return ave;
}

double MYMATH_sqave(double *data,int numstep) {
  int i;
  double sqave;
  for (i=0;i<numstep;++i) sqave+=data[i]*data[i];
  sqave=sqave/numstep;

  return sqave;
}

double MYMATH_var(double *data,int numstep) {
  double ave,sqave;
  double var;

  MYMATH_ave(data,numstep);
  MYMATH_sqave(data,numstep);

  return sqave-ave*ave;
}

double calc_flu(int num, double *data) {
  return sqrt(calc_var(num,data));
}

double max_f(int num,double *data) {
  int i;
  double max;
  max=data[0];
  for (i=1;i<num;++i)
    if (max < data[i]) max=data[i];
  return max;
}

double min_f(int num,double *data) {
  int i;
  double min;
  min=data[0];
  for (i=1;i<num;++i)
    if (min > data[i]) min=data[i];
  return min;
}

void mac_min_f(int num,double *data, double *max, double *min) {
  int i;
  *max=data[0];
  *min=*max;
  for (i=1;i<num;++i) {
    if (*max < data[i]) *max=data[i];
    if (*min > data[i]) *min=data[i];
  }
}

double ts_normalize(double *timeseries, int nums, int numv, double *timeseriesnorm) {
  int i,j;
  double *ave,*var;

  ave=(double *)gcemalloc(sizeof(double)*numv);
  var=(double *)gcemalloc(sizeof(double)*numv);

  for (i=0;i<numv;++i)
    ave[i]=0.0;
  for (i=0;i<numv;++i)
    for (j=0;j<nums;++j)
      ave[i]=(j*ave[i]+timeseries[j*numv+i])/(j+1);
  for (i=0;i<nums;++i)
    for (j=0;j<numv;++j)
      timeseriesnorm[i*numv+j]=timeseries[i*numv+j]-ave[j];
  for (i=0;i<numv;++i)
    var[i]=0.0;
  for (i=0;i<nums;++i)
    for (j=0;j<numv;++j)
      var[j]+=timeseriesnorm[i*numv+j]*timeseriesnorm[i*numv+j];

  for (j=0;j<numv;++j)
    var[j]=sqrt(var[j]/nums);
  for (i=0;i<nums;++i)
    for (j=0;j<numv;++j)
      timeseriesnorm[i*numv+j]=timeseriesnorm[i*numv+j]/var[j];

  for (i=0;i<numv;++i)
    var[i]=0.0;
  for (i=0;i<nums;++i)
    for (j=0;j<numv;++j)
      var[j]+=timeseriesnorm[i*numv+j]*timeseriesnorm[i*numv+j];

  for (j=0;j<numv;++j)
    var[j]=sqrt(var[j]/nums);

}

double *ts_normalize_od(double *timeseries,int nums) {
  int i,j;
  double ave,var;
  double *timeseriesnorm;

  timeseriesnorm=(double *)gcemalloc(sizeof(double)*nums);

  ave=0.0;
  for (i=0;i<nums;++i) ave=(i*ave+timeseries[i])/(i+1);
  for (i=0;i<nums;++i) timeseriesnorm[i]=timeseries[i]-ave;

  var=0.0;
  for (i=0;i<nums;++i) var+=timeseriesnorm[i]*timeseriesnorm[i];

  var=sqrt(var/nums);
  for (i=0;i<nums;++i) timeseriesnorm[i]=timeseriesnorm[i]/var;

  var=0.0;
  for (i=0;i<nums;++i) var+=timeseriesnorm[i]*timeseriesnorm[i];

  var=sqrt(var/nums);

  return timeseriesnorm;
}

double ts_normalize_v2(double **timeseries, int nums, int numv, double **timeseriesnorm) {
  int i,j;
  double *ave,*var;

  ave=(double *)gcemalloc(sizeof(double)*numv);
  var=(double *)gcemalloc(sizeof(double)*numv);

  for (i=0;i<numv;++i) ave[i]=0.0;
  for (i=0;i<numv;++i) for (j=0;j<nums;++j) ave[i]=(j*ave[i]+timeseries[i][j])/(j+1);
  for (i=0;i<numv;++i) for (j=0;j<nums;++j) timeseriesnorm[i][j]=timeseries[i][j]-ave[i];
  for (i=0;i<numv;++i) var[i]=0.0;
  for (i=0;i<numv;++i) for (j=0;j<nums;++j) var[i]+=timeseriesnorm[i][j]*timeseriesnorm[i][j];

  for (i=0;i<numv;++i) var[i]=sqrt(var[i]/nums);
  for (i=0;i<numv;++i) for (j=0;j<nums;++j) timeseriesnorm[i][j]=timeseriesnorm[i][j]/var[i];

  // for check
  for (i=0;i<numv;++i) var[i]=0.0;
  for (i=0;i<numv;++i) for (j=0;j<nums;++j)
      var[i]+=timeseriesnorm[i][j]*timeseriesnorm[i][j];

  for (i=0;i<numv;++i) var[i]=sqrt(var[i]/nums);
  //
}

int ts_bl_ave(double **timeseries, int nts, int numv, int nbts, int nvts, double **timeseriesba) {
  int i,j,k,numpoints;
  double *ave;
  
  numpoints=(int)nts/nvts;
  for (i=0;i<numv;++i) timeseriesba[i]=(double *)gcerealloc(timeseriesba[i],sizeof(double)*numpoints);
  ave=(double *)gcemalloc(sizeof(double)*numv);
  
  for (k=0;k<numv;++k) {
    for (i=0;i<numpoints;++i) {
      ave[k]=0.0;
      for (j=i*nvts;j<(i*nvts+nbts);++j)
	ave[k]=((j-i*nvts)*ave[k]+timeseries[k][j])/(j-i*nvts+1);
      timeseriesba[k][i]=ave[k];
    }
  }
  return numpoints;
}

double mm_len(double *vec) {
  int i;
  double len=0.0;

  for (i=0;i<3;++i) len+=vec[i]*vec[i];

  return sqrt(len);
}

int mm_outprod(double v1[3],double v2[3],double *v1x2) {

  v1x2[0]=v1[1]*v2[2]-v1[2]*v2[1];
  v1x2[1]=v1[2]*v2[0]-v1[0]*v2[2];
  v1x2[2]=v1[0]*v2[1]-v1[1]*v2[0];

}



/****************************************************************************************************************/
/* double ts_tra_cos_sin(double *timeseries, int nts, int numv, double *timeseries_sin_cos, int flag) {        */
/*   int i,j;												        */
/*   double pi;												        */
/* 													        */
/*   pi=acos(-1.0);											        */
/* 													        */
/*   for (i=0;i<nts;++i) {										        */
/*     for (j=0;j<numv;++j) {										        */
/*       if (flag=='R') {										        */
/* 	timeseries_sin_cos[i*numv*2+j*2]=cos(timeseries[i*numv+j]);					        */
/* 	timeseries_sin_cos[i*numv*2+j*2+1]=sin(timeseries[i*numv+j]);					        */
/*       }												        */
/*       else {												        */
/* 	timeseries_sin_cos[i*numv*2+j*2]=cos(timeseries[i*numv+j]/180.0*pi);				        */
/* 	timeseries_sin_cos[i*numv*2+j*2+1]=sin(timeseries[i*numv+j]/180.0*pi);				        */
/*       }												        */
/*     }												        */
/*   }													        */
/* }													        */
/****************************************************************************************************************/

