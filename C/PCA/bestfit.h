
#include "const.h"

//#define pi 3.14159265
//#define MAXNUMATOM 100

double bestfit(double coord_ref[MAXNUMATOM][3]
	       ,double coord_tag[MAXNUMATOM][3]
	       ,double velo_tag[MAXNUMATOM][3]
	       ,double mass[MAXNUMATOM],  int numatom);

void bestfit_iteration(double *traj,
		       double *crd_ave,
		       int numiteration,
		       int time,
		       int numatom,
		       double mass[MAXNUMATOM]);

