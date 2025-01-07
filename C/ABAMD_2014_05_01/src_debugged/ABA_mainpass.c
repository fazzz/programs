
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "ABA.h" // 2014-06-18
#include "ABAb.h"  // 2014-06-18
#include "LA.h"
#include "EF.h" // 2014-06-17

void ABAm_mainpass(ABI* abi,CLT *clt,double *Q,int numclut) {
  int i,j,k,nb;
  int nparent;
  double TMI[4][6][6],preABII[4][6][6],preBFI[4][6];

  for (i=numclut-1;i>=0;--i) {
    if (clt[i].terminal==TERM) {
       for (nb=0;nb<4;++nb) {
	for (j=0;j<6;++j) {
	  preBFI[nb][j]=0.0;
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=0.0;
	    preABII[nb][j][k]=0.0;
	  }
	}
      }
      ABAm_cCABI(abi[i].corABI,preABII,clt[i].IM,TMI,clt[i].num_branch);
    }
    else {
      for (nb=0;nb<clt[i].num_branch;++nb) {
	for (j=0;j<6;++j) {
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=clt[clt[i].nNumClutOfChild[nb]-1].TM[j][k];
	    preABII[nb][j][k]=abi[clt[i].nNumClutOfChild[nb]-1].preABI[j][k];
	  }
	}
      }
      ABAm_cCABI(abi[i].corABI,preABII,clt[i].IM,TMI,clt[i].num_branch);
    }

    abi[i].D=ABAm_cD(abi[i].corABI);
    ABAm_cKG(abi[i].KG,abi[i].corABI);
    ABAm_cPABI(abi[i].preABI,abi[i].corABI);

    if (clt[i].terminal==TERM) {
      for (nb=0;nb<4;++nb) {
	for (j=0;j<6;++j) {
	  preBFI[nb][j]=0.0;
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=0.0;
	  }
	}
      }
      //      for(nb=0;nb<clt[i].num_branch;++nb)  // 2014-08-13
      //	for (j=0;j<6;++j)  // 2014-08-13
      //	  cos(preBFI[nb][j]); // 2014-08-13
      ABAm_corBF(abi[i].corBF,preBFI,abi[i].corABI,clt[i].Coacc,clt[i].Cofrc,clt[i].Spfrc,TMI,clt[i].num_branch);
      //      for(nb=0;nb<clt[i].num_branch;++nb)  // 2014-08-13
      //	for (j=0;j<6;++j)  // 2014-08-13
      //	  cos(preBFI[nb][j]); // 2014-08-13
      //      for (j=0;j<6;++j) cos(abi[i].corBF[j]); // 2014-08-13
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
      //      for(nb=0;nb<clt[i].num_branch;++nb)  // 2014-08-13
      //	for (j=0;j<6;++j)  // 2014-08-13
      //	  cos(preBFI[nb][j]); // 2014-08-13
      ABAm_corBF(abi[i].corBF,preBFI,abi[i].corABI,clt[i].Coacc,clt[i].Cofrc,clt[i].Spfrc,TMI,clt[i].num_branch);
      //      for(nb=0;nb<clt[i].num_branch;++nb)  // 2014-08-13
      //	for (j=0;j<6;++j)  // 2014-08-13
      //	  cos(preBFI[nb][j]); // 2014-08-13
      //      for (j=0;j<6;++j) cos(abi[i].corBF[j]); // 2014-08-13
    }
    abi[i].eata=ABAm_eata(Q[i],abi[i].corBF);
    abi[i].nyu=ABAm_nyu(abi[i].eata,abi[i].D);
    ABAm_preBF(abi[i].preBF,abi[i].corBF,abi[i].KG,abi[i].eata);
  }
}

void ABAm_mainpass_TermOn(ABI* abi,CLT *clt,double *Q,int numclut,double *acc_Term,double *acc_Term2,double *vel_Term) {
  int i,j,k,nb;
  int nparent;
  double TMI[4][6][6],preABII[4][6][6],preBFI[4][6];

  for (i=numclut-1;i>=0;--i) {
    if (clt[i].terminal==TERM) {
       for (nb=0;nb<4;++nb) {
	for (j=0;j<6;++j) {
	  preBFI[nb][j]=0.0;
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=0.0;
	    preABII[nb][j][k]=0.0;
	  }
	}
      }
      ABAm_cCABI(abi[i].corABI,preABII,clt[i].IM,TMI,clt[i].num_branch);
    }
    else {
      for (nb=0;nb<clt[i].num_branch;++nb) {
	for (j=0;j<6;++j) {
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=clt[clt[i].nNumClutOfChild[nb]-1].TM[j][k];
	    preABII[nb][j][k]=abi[clt[i].nNumClutOfChild[nb]-1].preABI[j][k];
	  }
	}
      }
      ABAm_cCABI(abi[i].corABI,preABII,clt[i].IM,TMI,clt[i].num_branch);
    }

    abi[i].D=ABAm_cD(abi[i].corABI);
    ABAm_cKG(abi[i].KG,abi[i].corABI);
    ABAm_cPABI(abi[i].preABI,abi[i].corABI);

    if (clt[i].terminal==TERM) {
      for (nb=0;nb<4;++nb) {
	for (j=0;j<6;++j) {
	  preBFI[nb][j]=0.0;
	  for (k=0;k<6;++k) {
	    TMI[nb][j][k]=0.0;
	  }
	}
      }
      ABAm_corBF(abi[i].corBF,preBFI,abi[i].corABI,clt[i].Coacc,clt[i].Cofrc,clt[i].Spfrc,TMI,clt[i].num_branch);
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
      ABAm_corBF(abi[i].corBF,preBFI,abi[i].corABI,clt[i].Coacc,clt[i].Cofrc,clt[i].Spfrc,TMI,clt[i].num_branch);
    }
    abi[i].eata=ABAm_eata(Q[i],abi[i].corBF);
    abi[i].nyu=ABAm_nyu(abi[i].eata,abi[i].D);
    ABAm_preBF(abi[i].preBF,abi[i].corBF,abi[i].KG,abi[i].eata);
  }
  ABAm_corBF_TERM(abi[0].corABI,abi[0].corBF,acc_Term,acc_Term2,vel_Term);
}

int ABAm_cPABI(double preABI[6][6],double corABI[6][6]) {
  int i,j;

  for (i=0;i<6;++i)
    for (j=0;j<6;++j)
      preABI[i][j]=corABI[i][j]-corABI[i][2]*corABI[2][j]/corABI[2][2];

  return 0; // 2014-07-04
}

void ABAm_cKG(double KG[6],double corABI[6][6]){
  int i;
  for(i=0;i<6;++i) KG[i]= corABI[i][2]/corABI[2][2];
}

double ABAm_cD(double corABI[6][6]){
  return corABI[2][2];
}

void ABAm_cCABI(double corABI[6][6],double preABI[4][6][6],double IM[6][6], double TM[4][6][6],int numbranch){
  int i,j,k,l,nb;

  for(i=0;i<6;++i) for(j=0;j<6;++j) corABI[i][j]=0.0;

  for (nb=0;nb<numbranch;++nb)
    for(i=0;i<6;++i)
      for(j=0;j<6;++j)
	for(k=0;k<6;++k)
	  for(l=0;l<6;++l)
	    corABI[i][j]+= TM[nb][i][k]*preABI[nb][k][l]*TM[nb][j][l];
  for(i=0;i<6;++i) for(j=0;j<6;++j) corABI[i][j]+=IM[i][j];
}

void ABAm_preBF(double preBF[6],double corBF[6],double KG[6], double eata){
  int i;

  for(i=0;i<6;++i) preBF[i]=KG[i]*eata+corBF[i];
}

double ABAm_eata(double T,double Corbf[6]){

  return T - Corbf[2];
}

double ABAm_nyu(double eata,double D){

  return eata/D;
}

void ABAm_corBF(double corBF[6],double preBF[4][6],double corABI[6][6],double Coracc[6],double Corfrc[6],double frc[6],double TM[4][6][6],int numbranch){
   int i,j,nb;

   //  for (i=0;i<6;++i) cos(corBF[i]); // 2014-08-13
  for(i=0;i<6;++i) corBF[i] = 0.0;
  //  for (i=0;i<6;++i) cos(corBF[i]); // 2014-08-13
  for(nb=0;nb<numbranch;++nb)
    for(i=0;i<6;++i)
      for(j=0;j<6;++j)
	corBF[i] += TM[nb][i][j]*preBF[nb][j];
  //  for(nb=0;nb<numbranch;++nb) { // 2014-08-13
  //    for (i=0;i<6;++i) { // 2014-08-13
  //      cos(preBF[nb][i]); // 2014-08-13
  //      for(j=0;j<6;++j) { // 2014-08-13
  //	cos(TM[nb][i][j]); // 2014-08-13
  //      } // 2014-08-13
  //    } // 2014-08-13
  //  } // 2014-08-13
  //  for (i=0;i<6;++i) cos(corBF[i]); // 2014-08-13

  for(i=0;i<6;++i)
    for(j=0;j<6;++j)
      corBF[i]+=corABI[i][j]*Coracc[j];
  for(i=0;i<6;++i) corBF[i]+=Corfrc[i];
  for(i=0;i<6;++i) corBF[i]-=frc[i];
  //  for (i=0;i<6;++i) cos(corBF[i]); // 2014-08-13
}

void ABAm_corBF_TERM(double corABI[6][6],double corBF[6],double *acc_Term,double *acc_Term2,double *vel_Term){
  int i,j,k;
  double *temp,*inv,add[6];

  //  temp=(double *)gcemalloc(sizeof(double)*6*6); // 2014-07-18
  //  inv=(double *)gcemalloc(sizeof(double)*6*6);  // 2014-07-18
  temp=(double *)gcemalloc(sizeof(double)*6*6); // 2014-07-18
  inv=(double *)gcemalloc(sizeof(double)*6*6);  // 2014-07-18

  for (i=0;i<6;++i)
    for (j=0;j<6;++j)
      temp[i*6+j]=corABI[i][j];
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

    free(temp); // 2014-07-18
    free(inv);  // 2014-07-18
}


