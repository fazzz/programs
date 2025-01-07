
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABAb.h"
//#include "EF.h"

double ABAbs_forc(CLTb *clt,double *frc,int numclut) {
  int i,j,k,l;
  int natomtotal;
  double ***f_atom,**f_clust,***f_atom_trans, **f_clust_trans;

  int nNumClut;
  int na,numatomstart;
  int nNumClut_check;

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
  
  for(i=0;i<numclut;++i) for(j=0;j<6;++j) clt[i].Spfrc[j] = 0.0;
  
  for(i=0;i<numclut;++i){
    for(k=0;k<3;++k){
      f_clust[i][k]=0.0;
      f_clust_trans[i][k]=0.0;
    }
    for(j=0;j<clt[i].num_atom_clust;++j) {
      for(k=0;k<3;++k){
	f_atom[i][j][k]=0.0;
	f_atom_trans[i][j][k]=0.0;
      }
    }
  }

  /*********************************************************/
  /* natomtotal=0;					   */
  /* for(i=0;i<numclut;++i){				   */
  /*   for(j=0;j<clt[i].num_atom_clust;++j) {		   */
  /*     for(k=0;k<3;++k){				   */
  /* 	f_clust[i][k]+=frc[(j+natomtotal)*3+k];		   */
  /* 	f_atom[i][j][k]=frc[(j+natomtotal)*3+k];	   */
  /*     }						   */
  /*   }						   */
  /*   natomtotal+=clt[i].num_atom_clust;		   */
  /* }							   */
  /*********************************************************/

  for(nNumClut=0;nNumClut<numclut;++nNumClut){
    if (nNumClut==0)
      numatomstart=0;
    else
      numatomstart=clt[nNumClut].origin_atom_a;
    na=0;
    j=1;
    natomtotal=0;
    for(i=0;i<clt[nNumClut].num_atom_clust;++i) {
      if (numatomstart+i+na==clt[nNumClut+j].origin_atom_a) {
	//	na=0;
	nNumClut_check=nNumClut+j;
	for (;;) {
	  na+=clt[nNumClut+j].num_atom_clust;
	  if ((clt[nNumClut+j].join==clt[nNumClut_check].join) && (clt[nNumClut+j].terminal==0)) {
	    break;
	  }
	  ++j;
	}
	++j;
      } 
      if (nNumClut==0) {
	for(l=0;l<3;++l){                      
	  f_clust[nNumClut][l]+=frc[(i+na)*3+l]; 
	  f_atom[nNumClut][i][l]=frc[(i+na)*3+l];
	}                                      
      }
      else {
	for(l=0;l<3;++l){                      
	  f_clust[nNumClut][l]+=frc[(clt[nNumClut].origin_atom_a-1+i+na)*3+l]; 
	  f_atom[nNumClut][i][l]=frc[(clt[nNumClut].origin_atom_a-1+i+na)*3+l];
	}                                      
      }
    }
  }
  
  for(i=0;i<numclut;++i)
    for (j=0;j<3;++j) for (k=0;k<3;++k) clt[i].Spfrc[j+3] += clt[i].trans_A_to_CN[j][k]*f_clust[i][k];

  for(i=0;i<numclut;++i)
    for(j=0;j<clt[i].num_atom_clust;++j)
      for (k=0;k<3;++k)
	for (l=0;l<3;++l)
	  f_atom_trans[i][j][k] += clt[i].trans_A_to_CN[k][l]*f_atom[i][j][l];

  /**********************************************************************************************************************/
  /* for(i=0;i<numclut;++i) {											        */
  /*   for(j=0;j<clt[i].num_atom_clust;++j) {									        */
  /*     clt[i].Spfrc[0] += -clt[i].xoord[j*3+2]*f_atom_trans[i][j][1]+clt[i].xoord[j*3+1]*f_atom_trans[i][j][2];       */
  /*     													        */
  /*     clt[i].Spfrc[1] +=  clt[i].xoord[j*3+2]*f_atom_trans[i][j][0]-clt[i].xoord[j*3+0]*f_atom_trans[i][j][2];       */
  /*     													        */
  /*     clt[i].Spfrc[2] += -clt[i].xoord[j*3+1]*f_atom_trans[i][j][0]+clt[i].xoord[j*3+0]*f_atom_trans[i][j][1];       */
  /*   }													        */
  /* }														        */
  /**********************************************************************************************************************/

  for(nNumClut=0;nNumClut<numclut;++nNumClut){
    if (nNumClut==0)
      numatomstart=0;
    else
      numatomstart=clt[nNumClut].origin_atom_a;
    na=0;
    j=1;
    natomtotal=0;
    for(i=0;i<clt[nNumClut].num_atom_clust;++i) {
      if (numatomstart+i+na==clt[nNumClut+j].origin_atom_a) {
	nNumClut_check=nNumClut+j;
	for (;;) {
	  na+=clt[nNumClut+j].num_atom_clust;
	  if ((clt[nNumClut+j].join==clt[nNumClut_check].join) && (clt[nNumClut+j].terminal==0)) {
	    break;
	  }
	  ++j;
	}
	++j;
      } 
      clt[nNumClut].Spfrc[0] += -clt[nNumClut].xoord[(i+na)*3+2]*f_atom_trans[nNumClut][i][1]
                        	+clt[nNumClut].xoord[(i+na)*3+1]*f_atom_trans[nNumClut][i][2];
      
      clt[nNumClut].Spfrc[1] +=  clt[nNumClut].xoord[(i+na)*3+2]*f_atom_trans[nNumClut][i][0]
                        	-clt[nNumClut].xoord[(i+na)*3+0]*f_atom_trans[nNumClut][i][2];
      
      clt[nNumClut].Spfrc[2] += -clt[nNumClut].xoord[(i+na)*3+1]*f_atom_trans[nNumClut][i][0]
	                        +clt[nNumClut].xoord[(i+na)*3+0]*f_atom_trans[nNumClut][i][1]; 
    }
  }
}

