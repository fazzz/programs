
#include "LEASQDV.h"

double least_sqrt_devi(double *data, int numstep, double *a, double *b) {
  int i;
  double sum_x=0.0,sum_x_2=0.0,sum_xy=0.0,sum_y=0.0;

  for (i=0;i<numstep;++i) {
    sum_x   += i;
    sum_xy  += data[i]*i;
    sum_x_2 += i*i;
    sum_y   += data[i];
  }

  *a=1.0/(sum_x_2*numstep-sum_x*sum_x)*(numstep*sum_xy-sum_x*sum_y);
  *b=1.0/(sum_x_2*numstep-sum_x*sum_x)*(-sum_x*sum_xy+sum_x_2*sum_y);

  return *a;
}
