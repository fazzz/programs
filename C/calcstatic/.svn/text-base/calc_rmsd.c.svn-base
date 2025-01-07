#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double calc_rmsd(double *data, int num_initial, int num_final) {

  int i,numstep;
  double val_initial;
  double ave=0.0;
  double var=0.0;
  double rmsd;

  numstep=num_final-num_initial;
  val_initial = data[num_initial];
  for (i=num_initial;i<num_final;++i) {
    ave += data[i];
    var += data[i]*data[i];
  }

  ave = ave/numstep;
  var = var/numstep;
  var = var - ave*ave;

  if (rmsd != 0.0) {
    rmsd = sqrt(var);
  }
  else {
    rmsd = 0.0;
  }

  return rmsd;
}
