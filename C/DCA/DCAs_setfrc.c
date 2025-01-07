
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "DCA.h"
#include "EF.h"
#include "LA.h"

double DCAs_forc(CLT *clt,double *frc,int numclut) {
  int i,j,k,l;
  int natomtotal;
  double ***f_atom,**f_clust,***f_atom_trans, **f_clust_trans;

  f_clust=(double **)gcemalloc(sizeof(double *)*numclut);
  f_clust_trans=(double **)gcemalloc(sizeof(double *)*numclut);
  for (i=0;i<numclut;++i) {
    f_clust[i]=(double *)gcemalloc(sizeof(double)*3);
    f_clust_trans[i]=(double *)gcemalloc(sizeof(double)*3);
  }
  f_atom=(double ***)gcemalloc(sizeof(double **)*numclut);
  f_atom_trans=(double ***)gcemalloc(sizeof(double **)*numclut);
  for (i=0;i<numclut;++i) {
    f_atom[i]=(double **)gcemalloc(sizeof(double *)*clt[i].num_atom_clust);
    f_atom_trans[i]=(double **)gcemalloc(sizeof(double *)*clt[i].num_atom_clust);
    for (j=0;j<clt[i].num_atom_clust;++j) {
      f_atom[i][j]=(double *)gcemalloc(sizeof(double)*3);
      f_atom_trans[i][j]=(double *)gcemalloc(sizeof(double)*3);
    }
  }
  
  for(i=0;i<numclut;++i) 
    for(j=0;j<6;++j)
      clt[i].Spfrc[j] = 0.0;
  
  for(i=0;i<numclut;++i){
    for(j=0;j<clt[i].num_atom_clust;++j) {
      f_clust[i][k]=0.0;
      f_clust_trans[i][k]=0.0;
      for(k=0;k<3;++k){
	f_atom[i][j][k]=0.0;
	f_atom_trans[i][j][k]=0.0;
      }
    }
  }

  natomtotal=0;
  for(i=0;i<numclut;++i){
    for(j=0;j<clt[i].num_atom_clust;++j) {
      for(k=0;k<3;++k){
	f_clust[i][k]+=frc[(j+natomtotal)*3+k];
	f_atom[i][j][k]=frc[(j+natomtotal)*3+k];
      }
    }
    natomtotal+=clt[i].num_atom_clust;
  }
  
  for(i=0;i<numclut;++i)
    for (j=0;j<3;++j) 
      clt[i].Spfrc[j] = 0.0;
  for(i=0;i<numclut;++i)
    for (j=0;j<3;++j) 
      for (k=0;k<3;++k) 
	clt[i].Spfrc[j] += clt[i].trans_A_to_CN[j][k]*f_clust[i][k];

  for(i=0;i<numclut;++i)
    for(j=0;j<clt[i].num_atom_clust;++j)
      for (k=0;k<3;++k)
	for (l=0;l<3;++l)
	  f_atom_trans[i][j][k] += clt[i].trans_A_to_CN[j][k]*f_atom[i][j][l];

  for(i=0;i<numclut;++i) {
    for(j=0;j<clt[i].num_atom_clust;++j) {
      clt[i].Spfrc[3]
	+= -clt[i].xoord[j*3+2]*f_atom_trans[i][j][1]
	+clt[i].xoord[j*3+1]*f_atom_trans[i][j][2];
      
      clt[i].Spfrc[4]
	+= clt[i].xoord[j*3+2]*f_atom_trans[i][j][0]
	-clt[i].xoord[j*3+0]*f_atom_trans[i][j][2];
      
      clt[i].Spfrc[5]
	+= -clt[i].xoord[j*3+1]*f_atom_trans[i][j][0]
	+clt[i].xoord[j*3+0]*f_atom_trans[i][j][1];      
    }
  }
}

