
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABAb.h"

void ABAb_Inverse_backpass(ABIb* abi,CLTb *clt,double *Q,int numclut) {
  int i,j,k,nb;
  int nparent;
  double TMI[4][6][6],Fpart[4][6];

  for (i=numclut-1;i>=0;--i) {
    if (clt[i].terminal==TERM) {
       for (nb=0;nb<4;++nb) {
	for (j=0;j<6;++j) {
	  Fpart[nb][j]=0.0;
	  for (k=0;k<6;++k)
	    TMI[nb][j][k]=0.0;
	}
      }
      Q[i]=ABAbb_cQ(clt[i].F,clt[i].IM,clt[i].Spacc,clt[i].Cofrc,clt[i].Spfrc,TMI,Fpart,clt[i].num_branch);
    }
    else {
      for (nb=0;nb<clt[i].num_branch;++nb) {
	for (j=0;j<6;++j) {
	  Fpart[nb][j]=clt[clt[i].nNumClutOfChild[nb]-1].F[j];
	  for (k=0;k<6;++k)
	    TMI[nb][j][k]=clt[clt[i].nNumClutOfChild[nb]-1].TM[j][k];
	}
      }
      Q[i]=ABAbb_cQ(clt[i].F,clt[i].IM,clt[i].Spacc,clt[i].Cofrc,clt[i].Spfrc,TMI,Fpart,clt[i].num_branch);
    }
  }
}

double ABAbb_cQ(double F[6], double IM[6][6], double Spacc[6],double Cofrc[6],double Spfrc[6], double TM[4][6][6],double Fparent[4][6],int num_branch){
  int i,j,nb;
  double Q;

  for(i=0;i<6;++i) F[i] = 0.0;
  for(nb=0;nb<num_branch;++nb)
    for(i=0;i<6;++i)
      for(j=0;j<6;++j)
	F[i] += TM[nb][i][j]*Fparent[nb][j];

  for(i=0;i<6;++i)
    for(j=0;j<6;++j)
      F[i]+=IM[i][j]*Spacc[j];
  for(i=0;i<6;++i) F[i]+=Cofrc[i];
  for(i=0;i<6;++i) F[i]-=Spfrc[i];

  Q=F[2];
  return Q;
}

