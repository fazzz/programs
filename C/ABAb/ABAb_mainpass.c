
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABAb.h"
#include "LA.h"

void ABAbm_mainpass(ABIb* abi,CLTb *clt,double *Q,int numclut) {
  int i,j,k,nb;
  int nparent;
  double TMI[10][6][6],preABIbI[10][6][6],preBFI[10][6];

  for (i=numclut-1;i>=0;--i) {
    if (clt[i].terminal==TERM) {
       for (nb=0;nb<10;++nb) {
	for (j=0;j<6;++j) {
	  preBFI[nb][j]=0.0;
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=0.0;
	    preABIbI[nb][j][k]=0.0;
	  }
	}
      }
      ABAbm_cCABIb(abi[i].corABIb,preABIbI,clt[i].IM,TMI,clt[i].num_branch);
    }
    else {
      for (nb=0;nb<clt[i].num_branch;++nb) {
	for (j=0;j<6;++j) {
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=clt[clt[i].nNumClutOfChild[nb]-1].TM[j][k];
	    preABIbI[nb][j][k]=abi[clt[i].nNumClutOfChild[nb]-1].preABIb[j][k];
	  }
	}
      }
      ABAbm_cCABIb(abi[i].corABIb,preABIbI,clt[i].IM,TMI,clt[i].num_branch);
    }

    abi[i].D=ABAbm_cD(abi[i].corABIb);
    ABAbm_cKG(abi[i].KG,abi[i].corABIb);
    ABAbm_cPABIb(abi[i].preABIb,abi[i].corABIb);

    if (clt[i].terminal==TERM) {
      for (nb=0;nb<10;++nb) {
	for (j=0;j<6;++j) {
	  preBFI[nb][j]=0.0;
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=0.0;
	  }
	}
      }
      ABAbm_corBF(abi[i].corBF,preBFI,abi[i].corABIb,clt[i].Coacc,clt[i].Cofrc,clt[i].Spfrc,TMI,clt[i].num_branch);
    }
    else {
      for (nb=0;nb<clt[i].num_branch;++nb) {
	for (j=0;j<6;++j) {
	  preBFI[nb][j]=abi[clt[i].nNumClutOfChild[nb]-1].preBF[j];
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=clt[clt[i].nNumClutOfChild[nb]-1].TM[j][k];
	  }
	}
      }
      ABAbm_corBF(abi[i].corBF,preBFI,abi[i].corABIb,clt[i].Coacc,clt[i].Cofrc,clt[i].Spfrc,TMI,clt[i].num_branch);
    }
    abi[i].eata=ABAbm_eata(Q[i],abi[i].corBF);
    abi[i].nyu=ABAbm_nyu(abi[i].eata,abi[i].D);
    ABAbm_preBF(abi[i].preBF,abi[i].corBF,abi[i].KG,abi[i].eata);
  }
}

void ABAbm_mainpass_TermOn(ABIb* abi,CLTb *clt,double *Q,int numclut,double *acc_Term,double *acc_Term2,double *vel_Term) {
  int i,j,k,nb;
  int nparent;
  double TMI[10][6][6],preABIbI[10][6][6],preBFI[10][6];

  for (i=numclut-1;i>=0;--i) {
    if (clt[i].terminal==TERM) {
       for (nb=0;nb<10;++nb) {
	for (j=0;j<6;++j) {
	  preBFI[nb][j]=0.0;
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=0.0;
	    preABIbI[nb][j][k]=0.0;
	  }
	}
      }
      ABAbm_cCABIb(abi[i].corABIb,preABIbI,clt[i].IM,TMI,clt[i].num_branch);
    }
    else {
      for (nb=0;nb<clt[i].num_branch;++nb) {
	for (j=0;j<6;++j) {
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=clt[clt[i].nNumClutOfChild[nb]-1].TM[j][k];
	    preABIbI[nb][j][k]=abi[clt[i].nNumClutOfChild[nb]-1].preABIb[j][k];
	  }
	}
      }
      ABAbm_cCABIb(abi[i].corABIb,preABIbI,clt[i].IM,TMI,clt[i].num_branch);
    }

    abi[i].D=ABAbm_cD(abi[i].corABIb);
    ABAbm_cKG(abi[i].KG,abi[i].corABIb);
    ABAbm_cPABIb(abi[i].preABIb,abi[i].corABIb);

    if (clt[i].terminal==TERM) {
      for (nb=0;nb<10;++nb) {
	for (j=0;j<6;++j) {
	  preBFI[nb][j]=0.0;
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=0.0;
	  }
	}
      }
      ABAbm_corBF(abi[i].corBF,preBFI,abi[i].corABIb,clt[i].Coacc,clt[i].Cofrc,clt[i].Spfrc,TMI,clt[i].num_branch);
    }
    else {
      for (nb=0;nb<clt[i].num_branch;++nb) {
	for (j=0;j<6;++j) {
	  preBFI[nb][j]=abi[clt[i].nNumClutOfChild[nb]-1].preBF[j];
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=clt[clt[i].nNumClutOfChild[nb]-1].TM[j][k];
	  }
	}
      }
      ABAbm_corBF(abi[i].corBF,preBFI,abi[i].corABIb,clt[i].Coacc,clt[i].Cofrc,clt[i].Spfrc,TMI,clt[i].num_branch);
    }
    abi[i].eata=ABAbm_eata(Q[i],abi[i].corBF);
    abi[i].nyu=ABAbm_nyu(abi[i].eata,abi[i].D);
    ABAbm_preBF(abi[i].preBF,abi[i].corBF,abi[i].KG,abi[i].eata);
  }
  ABAbm_corBF_TERM(abi[0].corABIb,abi[0].corBF,acc_Term,acc_Term2,vel_Term);
}

int ABAbm_cPABIb(double preABIb[6][6],double corABIb[6][6]) {
  int i,j;

  for (i=0;i<6;++i)
    for (j=0;j<6;++j)
      preABIb[i][j]=corABIb[i][j]-corABIb[i][2]*corABIb[2][j]/corABIb[2][2];
}

void ABAbm_cKG(double KG[6],double corABIb[6][6]){
  int i;
  for(i=0;i<6;++i) KG[i]= corABIb[i][2]/corABIb[2][2];
}

double ABAbm_cD(double corABIb[6][6]){
  return corABIb[2][2];
}

void ABAbm_cCABIb(double corABIb[6][6],double ***preABIb/*[4][6][6]*/,double IM[6][6], double ***TM/*[4][6][6]*/,int numbranch){
  int i,j,k,l,nb;

  for(i=0;i<6;++i) for(j=0;j<6;++j) corABIb[i][j]=0.0;

  for (nb=0;nb<numbranch;++nb)
    for(i=0;i<6;++i)
      for(j=0;j<6;++j)
	for(k=0;k<6;++k)
	  for(l=0;l<6;++l)
	    corABIb[i][j]+= TM[nb][i][k]*preABIb[nb][k][l]*TM[nb][j][l];
  for(i=0;i<6;++i) for(j=0;j<6;++j) corABIb[i][j]+=IM[i][j];
}

/******************************************************************************************************************************************/
/* void ABAbm_cCABIb(double corABIb[6][6],double preABIb[/\*10*\/4][6][6],double IM[6][6], double TM[/\*10*\/4][6][6],int numbranch){	  */
/*   int i,j,k,l,nb;															  */
/* 																	  */
/*   for(i=0;i<6;++i) for(j=0;j<6;++j) corABIb[i][j]=0.0;										  */
/* 																	  */
/*   for (nb=0;nb<numbranch;++nb)													  */
/*     for(i=0;i<6;++i)															  */
/*       for(j=0;j<6;++j)														  */
/* 	for(k=0;k<6;++k)														  */
/* 	  for(l=0;l<6;++l)														  */
/* 	    corABIb[i][j]+= TM[nb][i][k]*preABIb[nb][k][l]*TM[nb][j][l];									  */
/*   for(i=0;i<6;++i) for(j=0;j<6;++j) corABIb[i][j]+=IM[i][j];										  */
/* }																	  */
/******************************************************************************************************************************************/

void ABAbm_preBF(double preBF[6],double corBF[6],double KG[6], double eata){
  int i;

  for(i=0;i<6;++i) preBF[i]=KG[i]*eata+corBF[i];
}

double ABAbm_eata(double T,double Corbf[6]){

  return T - Corbf[2];
}

double ABAbm_nyu(double eata,double D){

  return eata/D;
}

void ABAbm_corBF(double corBF[6],double **preBF/*[10][6]*/,double corABIb[6][6],double Coracc[6],double Corfrc[6],double frc[6],double ***TM/*[10][6][6]*/,int numbranch){
   int i,j,nb;

  for(i=0;i<6;++i) corBF[i] = 0.0;
  for(nb=0;nb<numbranch;++nb)
    for(i=0;i<6;++i)
      for(j=0;j<6;++j)
	corBF[i] += TM[nb][i][j]*preBF[nb][j];

  for(i=0;i<6;++i)
    for(j=0;j<6;++j)
      corBF[i]+=corABIb[i][j]*Coracc[j];
  for(i=0;i<6;++i) corBF[i]+=Corfrc[i];
  for(i=0;i<6;++i) corBF[i]-=frc[i];
}

void ABAbm_corBF_TERM(double corABIb[6][6],double corBF[6],double *acc_Term,double *acc_Term2,double *vel_Term){
  int i,j,k;
  double *temp,*inv,add[6];

  temp=(double *)gcemalloc(sizeof(double)*6*6);
  inv=(double *)gcemalloc(sizeof(double)*6*6);

  for (i=0;i<6;++i)
    for (j=0;j<6;++j)
      temp[i*6+j]=corABIb[i][j];
  invm2(temp,inv,6);
    for(i=0;i<6;++i) {
      acc_Term[i]=0.0;
      acc_Term2[i]=0.0;
    }
    add[0]=0.0;
    add[1]=0.0;
    add[2]=0.0;
    add[3]=(-vel_Term[2]*vel_Term[3+1]+vel_Term[1]*vel_Term[3+2]);
    add[4]=( vel_Term[2]*vel_Term[3+0]-vel_Term[0]*vel_Term[3+2]);
    add[5]=(-vel_Term[1]*vel_Term[3+0]+vel_Term[0]*vel_Term[3+1]);

    for(i=0;i<6;++i)
      for(j=0;j<6;++j)
	acc_Term2[i]-=inv[i*6+j]*corBF[j];
    for(i=0;i<6;++i) acc_Term[i]=acc_Term2[i]-add[i];
}


