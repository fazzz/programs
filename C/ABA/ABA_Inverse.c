
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABA.h"
#include "EF.h"

double solverABA_Inverse(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom){
  ABI *abi;

  abi=(ABI *)gcemalloc(sizeof(ABI)*numclt);

  ABAs_forc(clt,frc,numclt);

  ABA_Inverse_mainpass(clt,qacc,qvel,numclt,numatom,crd);
  ABA_Inverse_backpass(abi,clt,Q,numclt);

}
