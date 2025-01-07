
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABAb.h"
#include "EF.h"

double solverABAb_Inverse(double *qacc,double *qvel,CLTb *clt,double *Q,double *frc,double *crd,int numclt,int numatom){
  ABIb *abi;

  abi=(ABIb *)gcemalloc(sizeof(ABIb)*numclt);

  ABAbs_forc(clt,frc,numclt);

  ABAb_Inverse_mainpass(clt,qacc,qvel,numclt,numatom,crd);
  ABAb_Inverse_backpass(abi,clt,Q,numclt);

}
