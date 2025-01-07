#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double calc_least_sqrt_devi(double *xdata, double *ydata, int numdata) {
  int i,j;
  double sum_x=0.0,sum_x_2=0.0,sum_xy=0.0,sum_y=0.0;
  double a,b;

  for (i=0;i<numdata;++i) {
    sum_x   += xdata[i];
    sum_xy  += xdata[i]*ydata[i];
    sum_x_2 += xdata[i]*xdata[i];
    sum_y   += ydata[i];
  }

  a=1.0/(sum_x_2*numdata-sum_x*sum_x)*(numdata*sum_xy-sum_x*sum_y);
  b=1.0/(sum_x_2*numdata-sum_x*sum_x)*(-sum_x*sum_xy-sum_x_2*sum_y);

  return a;
}
