
#ifndef INCLUDE_LogMFD_cMF
#define INCLUDE_LogMFD_cMF

#include "FFLc.h"

double LogMFD_cdF_dX_wrtdihed(double *crd,struct potential *e, struct force *f, struct AmberParmL ap,
			      int numatom, int na1, int na2);

#endif
