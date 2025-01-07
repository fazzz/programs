
#include "DHAD.h"
#include "TOPO.h"

double dhad(int numd, int *list, double *crd, double *crdref) {
  int i,j,k;
  double pi,atom[4][3],atomref[4][3],diha,diharef,dhad,dhadtotal=0.0;

  pi=acos(-1.0);

  for (i=0;i<numd;++i) {
    for (j=0;j<4;++j) {
      for (k=0;k<3;++k) {
	atom[j][k]=crd[list[i*4+j]*3+k];
	atomref[j][k]=crdref[list[i*4+j]*3+k];
      }
    }

    diha=dih(atom[0],atom[1],atom[2],atom[3]);
    diharef=dih(atomref[0],atomref[1],atomref[2],atomref[3]);

    dhad=fabs(diha-diharef);
    if (dhad > 2.0*pi-dhad) dhad=2.0*pi-dhad;
    dhadtotal=(i*dhadtotal+dhad)/(i+1);
  }
  return dhadtotal/pi;

}

