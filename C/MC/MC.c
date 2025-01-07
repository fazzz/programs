
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "MC.h"
#include "RAND.h"
#include "EF.h"

int Metropolis(double delta) {
  double expdelta; // exp(-delta)
  double u;

  expdelta=exp(-1.0*delta);
  u=genrand_real2();

  if (expdelta >= 1.0)
    return 1;
  else if (u <= expdelta)
    return 1;
  else
    return 0;
}

