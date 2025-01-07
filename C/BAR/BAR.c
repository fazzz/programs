
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "BAR.h"
#include "mymath.h"

int BAR_ite(int n1,         // num of snapshots of U1
	    int n2,         // num of snapshots of U2
	    double C,       // C=b-1ln(Q0n1/Q1n0)
	    double **ebU1,  // exp(-bU1(Ri,j))
	    double **ebU2,  // exp(-bU2(Ri,j))
	    double dF       // d Free Enrgy
	    ) {
  int i,n;
  double comv;
  double C_new;
  double deno,nume;
  double delta;

  C=0.0;
  for (n=0;n<maxi;++n)  {
    deno=0.0;
    nume=0.0;
    for (i=0;i<n1;++i)
      deno+=Fermi(ebU1[0][i]-ebU2[0][i]+C);
    for (i=0;i<n1;++i)
      nume+=Fermi(ebU1[1][i]-ebU2[1][i]+C);
    dF=beta*log(deno/nume)+C;
    C_new=1.0/beta*log(n1/n0)-dF;
    
    conv=1;
    delta=abs(C_new-C);
    if (delta>=MINSHIFT) conv=0;
    ebF[k]=ebmF[0]/ebmF[k];
    if (conv>0) break;
  }

  return n;
}

int NBB_ite(int n1,         // num of snapshots of U1
	    int n2,         // num of snapshots of U2
	    double C,       // C=b-1ln(Q0n1/Q1n0)
	    double **ebU1,  // exp(-bU1(Ri,j))
	    double **ebU2,  // exp(-bU2(Ri,j))
	    double ***ebW,
	    double dF       // d Free Enrgy
	    ) {
  int i,n;
  double comv;
  double C_new;
  double deno,nume;
  double delta;

  C=1.0;
  for (n=0;n<maxi;++n)  {
    deno=0.0;
    nume=0.0;
    for (i=0;i<n1;++i)
      deno+=Fermi(ebU1[0][i]-ebU2[0][i]+C);
    for (i=0;i<n1;++i)
      nume+=Fermi(ebU1[1][i]-ebU2[1][i]+C);
    dF=beta*log(deno/nume)+C;
    C_new=1.0/beta*log(n1/n0)-dF;
    
    conv=1;
    delta=abs(C_new-C);
    if (delta>=MINSHIFT) conv=0;
    ebF[k]=ebmF[0]/ebmF[k];
    if (conv>0) break;
  }

  return n;
}
