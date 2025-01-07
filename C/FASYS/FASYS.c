#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "FASYS.h"
#include "EF.h"
#include "TOPO.h"

double FASYS_calcp(double crd[5][3],double *p,double n[2]) {
  int i;
  double pi;
  double v=0.0;
  double kb=100.0,ka=50.0;
  double kd1=0.0,n1=0.0,kd2=0.0,n2=0.0;
  double l_eq=1.53,a_eq;
  pi=acos(-1.0);
  a_eq=111.0/180*pi;

  (*energy).bond=0.0;
  (*energy).angl=0.0;
    
  for (i=0;i<4;++i)
    (*energy).bond+=0.5*kb*(len(crd[i],crd[i+1])-l_eq)*(len(crd[i],crd[i+1])-l_eq);

  for (i=0;i<3;++i)
    (*energy).angl+=0.5*ka*(ang(crd[i],crd[i+1],crd[i+2])-a_eq)*(ang(crd[i],crd[i+1],crd[i+2])-a_eq);

  kd1=2.5;
  n1=3.0;
  kd2=2.5;
  n2=3.0;
    /* kd1=2.0; */
    /***********/
    /* n1=3.0; */
    /***********/
    /* kd2=2.0; */
    /***********/
    /* n2=1.0; */
    /***********/

  *p1=dih(crd[0],crd[1],crd[2],crd[3]);
  *h1=dih(crd[1],crd[2],crd[3],crd[4]);

  (*energy).dihd=0.5*kd1*(1.0+cos(n1*(*p1)))+0.5*kd2*(1.0+cos(n2*(*h1)));
  (*energy).v=(*energy).bond+(*energy).angl+(*energy).dihd;

  v=(*energy).bond[0]+(*energy).angl[0]+(*energy).dihd[0]+(*energy).eata[0];

  return v;
}
