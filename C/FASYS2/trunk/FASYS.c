
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "FASYS.h"
#include "TOPO.h"
#include "QUA.h"

double set_coord_from_dihed(double phi, double psi, double q[5][3]) {
  int i,j,k;
  double pi;
  double u[3];
  double q_dummy[4],rotation[4],q_rotated[4];
  double l_eq=1.53,ang_eq;
  pi=acos(-1.0);
  ang_eq=111.0/180*pi;

  for (i=0;i<3;++i) q[0][i]=0.0;
  q[1][0]=l_eq;
  for (i=1;i<3;++i) q[1][i]=0.0;
  q[2][0]=cos(-ang_eq)*l_eq+q[1][0];
  q[2][1]=sin(-ang_eq)*l_eq+q[1][1];
  q[2][2]=0.0;
  q[3][0]=cos(-2.0*ang_eq)*l_eq+q[2][0];
  q[3][1]=sin(-2.0*ang_eq)*l_eq+q[2][1];
  q[3][2]=0.0;
  q[4][0]=cos(-3.0*ang_eq)*l_eq+q[3][0];
  q[4][1]=sin(-3.0*ang_eq)*l_eq+q[3][1];
  q[4][2]=0.0;

  q_dummy[0]=0.0;
  for (i=1;i<4;++i) q_dummy[i]=q[3][i-1]-q[1][i-1];
  for (i=0;i<3;++i) u[i]=(q[2][i]-q[1][i])/l_eq;
  rotation[0]=cos(0.5*phi);
  for (i=1;i<4;++i) rotation[i]=sin(0.5*phi)*u[i-1];
  qua_rot(q_dummy,rotation,q_rotated);
  for (i=0;i<3;++i) q[3][i]=q_rotated[i+1]+q[1][i];

  q_dummy[0]=0.0;
  for (i=1;i<4;++i) q_dummy[i]=q[4][i-1]-q[1][i-1];
  qua_rot(q_dummy,rotation,q_rotated);
  for (i=0;i<3;++i) q[4][i]=q_rotated[i+1]+q[1][i];

  q_dummy[0]=0.0;
  for (i=1;i<4;++i) q_dummy[i]=q[4][i-1]-q[2][i-1];
  for (i=0;i<3;++i) u[i]=(q[3][i]-q[2][i])/l_eq;
  rotation[0]=cos(0.5*psi);
  for (i=1;i<4;++i) rotation[i]=sin(0.5*psi)*u[i-1];
  qua_rot(q_dummy,rotation,q_rotated);
  for (i=0;i<3;++i) q[4][i]=q_rotated[i+1]+q[2][i];

}

double FASYS_calcpote(double q[5][3], double kd[2],double n[2]) {
  int i;
  double pi;
  double v=0.0;
  double kb=100.0,ka=50.0;
  double l_eq=1.53,a_eq;
  pi=acos(-1.0);
  a_eq=111.0/180*pi;

  for (i=0;i<4;++i) v+=0.5*kb*(len(q[i],q[i+1])-l_eq)*(len(q[i],q[i+1])-l_eq);
  for (i=0;i<3;++i) v+=0.5*ka*(ang(q[i],q[i+1],q[i+2])-a_eq)*(ang(q[i],q[i+1],q[i+2])-a_eq);
  
  v+=0.5*kd[0]*(1.0+cos(n[0]*(dih(q[0],q[1],q[2],q[3]))))+0.5*kd[1]*(1.0+cos(n[1]*dih(q[1],q[2],q[3],q[4])));

  return v;
}
