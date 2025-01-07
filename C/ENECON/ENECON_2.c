
#include "ENECON_2.h"

double ec_flu_2(double *ene, int numstep) {
  int i,j;
  double ave=0.0;
  double var=0.0;
  double val_i;
  double rmsd;

  double sqrt_var,fabs_val_i;

  val_i = ene[0];
  for (i=0;i<numstep;++i) {
    ave += ene[i];
  }
  ave = ave/numstep;

  for (i=0;i<numstep;++i) {
    var+=(ave - ene[i])*(ave - ene[i]);
  }
  var = var/numstep;

  if (var != 0.0 && val_i != 0.0) {
    sqrt_var=sqrt(var);
    fabs_val_i=fabs(val_i);
    rmsd = sqrt_var/fabs_val_i;
  }
  else if ( var==0.0 )
    rmsd = 0.0;
  else {
    printf("error\n");
    exit(1);
  }

  return rmsd;
}

double ec_avd_2(double *ene, int numstep) {
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



