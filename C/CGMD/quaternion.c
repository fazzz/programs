
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

void qua_trans_omgtodqua(double omg[3],double q[4],double dq[4]) {
  int i,j;
  double mat[4][4];
  double omgd[4],dq_temp[4];

  mat[0][0]=-q[3];
  mat[0][1]=-q[0];
  mat[0][2]= q[1];
  mat[0][3]= q[2];

  mat[1][0]= q[0];
  mat[1][1]=-q[3];
  mat[1][2]=-q[2];
  mat[1][3]= q[1];

  mat[2][0]= q[2];
  mat[2][1]= q[1];
  mat[2][2]= q[0];
  mat[2][3]= q[3];

  mat[3][0]=-q[1];
  mat[3][1]= q[2];
  mat[3][2]=-q[3];
  mat[3][3]= q[0];

  for (i=0;i<3;++i)
    omgd[i]=omg[i];
  omgd[3]=0.0;

  for (i=0;i<4;++i)
    dq_temp[i]=0.0;
  for (i=0;i<4;++i)
    for (j=0;j<4;++j)
      dq_temp[i]+=0.5*mat[i][j]*omgd[j];

  dq[0]=dq_temp[3];
  dq[1]=dq_temp[1];
  dq[2]=dq_temp[0];
  dq[3]=dq_temp[2];

}

void qua_set_rotmat(double q[4],double rotmat[3][3]) {

  rotmat[0][0] = q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
  rotmat[1][1] = q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
  rotmat[2][2] = q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
  
  rotmat[0][1] = 2.0*(q[1]*q[2]-q[0]*q[3]);
  rotmat[1][0] = 2.0*(q[1]*q[1]+q[0]*q[3]);
  
  rotmat[0][2] = 2.0*(q[1]*q[3]+q[0]*q[2]);
  rotmat[2][0] = 2.0*(q[1]*q[3]-q[0]*q[2]);
  
  rotmat[1][2] = 2.0*(q[2]*q[3]-q[0]*q[1]);
  rotmat[2][1] = 2.0*(q[2]*q[3]+q[0]*q[1]);

}

void qua_euler(double euler[3], double q[4]) {
  q[0]=cos(0.5*euler[1])*cos(0.5*(euler[0]+euler[1]));
  q[1]=sin(0.5*euler[1])*cos(0.5*(euler[0]-euler[1]));
  q[2]=sin(0.5*euler[1])*sin(0.5*(euler[0]-euler[1]));
  q[3]=cos(0.5*euler[1])*sin(0.5*(euler[0]+euler[1]));
}
