#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Rand.h"
#include "Metropolis.h"

int Metropolis(double Etrial, double E) {
  int i;
  double randnum;

  if (Etrial/E < 1.0) {
    return ACP;
  }
  else if (((randnum = randMT()) < Etrial/E) < 1.0) {
    return ACP;
  }
  else {
    return REJ;
  }
}
