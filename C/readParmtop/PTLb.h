#ifndef INCLUDE_PTLb
#define INCLUDE_PTLb

#include "PTL.h"

int readParmtopLb(FILE *parmfile,struct AmberParmL*ap);
int readdihedpairsLb(int **atomdihedpairs, int *num,struct AmberParmL ap);
int PTL_joinatomtoresb(int numatom, char LABERES[4],struct AmberParmL ap);

int PTL_res_cab(int numres,struct AmberParmL ap);

//int PTL_iniatomnumofres(int numres);
int PTL_resnumb(int numatom,struct AmberParmL ap);

int PTL_resnum2b(int numatom,struct AmberParmL ap);

int PTL_canum_fromresnumb(int numres,struct AmberParmL ap);

int PTL_which_includeb(int numres,int *listres,int numlistres,struct AmberParmL ap);

#endif
