
#include <stdio.h>
#include <math.h>
#include "quaternion.h"

void quaternion_rotation(double q[4],double r[4],double roted[4])
{
  int i,j;
  double qs[4],t[4];

  for (j=0;j<4;++j)
    roted[j] = 0.0;
  quaternion_conjugate(q, qs);
  quaternion_product(r, qs, t);
  quaternion_product(q, t, roted);
      
}

void quaternion_product(double r[4],double q[4], double p[4])
{
  int i,j;
  double tempmat[4][4];
  double out[3];

  p[0]=0.0;
  for (i=1;i<4;++i)
    p[0] -= r[i]*q[i];

  p[0]+=r[0]*q[0];

  out[0] = r[2]*q[3]-r[3]*q[2];
  out[1] = r[3]*q[1]-r[1]*q[3];
  out[2] = r[1]*q[2]-r[2]*q[1];

  for (i=1;i<4;++i)
    p[i] = q[0]*r[i]+r[0]*q[i]+out[i-1];


}

void quaternion_conjugate(double q[4],double qs[4])
{
   qs[0] =  q[0];
   qs[1] = -q[1];
   qs[2] = -q[2];
   qs[3] = -q[3];

}

