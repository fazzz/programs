
#ifndef INCLUDE_TA_MDb
#define INCLUDE_TA_MDb

#include <netcdf.h>
#include "netcdf_mineL.h"

#include "FFLc.h"

double ffLc_calcffandforce_AA(double *crd, int numatom,struct potential *ene,struct force *f,struct AmberParmL ap);

#endif
