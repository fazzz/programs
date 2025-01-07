
#include "ENECON.h"

double ec_flu(double *ene, int numstep) {
  int i,j;
  double ave=0.0;
  double var=0.0;
  double val_i;
  double rmsd;

  val_i = ene[0];
  for (i=0;i<numstep;++i) {
    ave += ene[i];
    var += ene[i]*ene[i];
  }

  ave = ave/numstep;
  var = var/numstep;
  var = var - ave*ave;

  if (var != 0.0 && val_i != 0.0)
    rmsd = sqrt(var)/fabs(val_i);
  else if ( var==0.0 )
    rmsd = 0.0;
  else {
    printf("error\n");
    exit(1);
  }

  return rmsd;
}

double ec_avd(double *ene, int numstep) {
  int i;
  double ave=0.0;

  if (ene[0] == 0.0) {
    printf("error\n");
    exit(1);
  }

  for (i=0;i<numstep;++i) ave += fabs(ene[i]-ene[0])/fabs(ene[0]);

  ave = ave/numstep;

  return ave;
}



