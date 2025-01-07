
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABAb.h"
#include "EF.h"
#include "LA.h"

void ABAb_Inverse_mainpass(CLTb *clt,double *qacc,double *qvel,int numclt,int numatom,double *crd) {
  int i,j;
  int nNumClutOfParent;

  for (i=0;i<numclt;++i) {
    nNumClutOfParent = clt[i].nNumClutOfParent-1;
    if (i==0) {
      for(j=0;j<6;++j) {
	clt[i].Spacc[j] = 0.0;
	clt[i].Spvel[j] = 0.0;
	clt[i].Cofrc[j] = 0.0;
	clt[i].Coacc[j] = 0.0;
      }
    }
    else {
      ABAbp_Spvelo(clt[i].Spvel,qvel[i],clt[i].TM,clt[nNumClutOfParent].Spvel);
    }
  }

  for (i=0;i<numclt;++i) {
    nNumClutOfParent = clt[i].nNumClutOfParent-1;
    ABAbp_Cofr(clt[i].Cofrc,clt[i].Spvel,clt[i].IM,clt[i].qCOM,clt[i].sum_mass);
    ABAb_Coacc(clt[i].Coacc,clt[i].Spvel,clt[nNumClutOfParent].Spvel,qvel[i],
	      clt[i].origin_atom_a,clt[nNumClutOfParent].origin_atom_a,nNumClutOfParent,
	      crd,clt[nNumClutOfParent].trans_A_to_CN,clt[i].TM);
    ABAb_Inverse_Spacc(clt[i].Spacc,qacc[i],clt[i].TM,clt[i].Coacc,clt[nNumClutOfParent].Spacc);
  }
}

void ABAb_Inverse_Spacc(double Spacc[6],double qacc,double TMat[6][6],double Coacc[6],double SpaccPart[6]) {
  int i,j;

  for(i=0;i<6;++i) Spacc[i] = 0.0;

  for(i=0;i<6;++i)
    for(j=0;j<6;++j)
      Spacc[i]+=TMat[j][i]*SpaccPart[j];

  for(i=0;i<6;++i) Spacc[i]+=Coacc[i];
  Spacc[2]+=qacc;
}

