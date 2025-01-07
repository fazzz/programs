#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "MD.h"
#include "BD.h"

double MatrixOfIJEle[4][4];

// 巨大行列の計算
void calc_Matrix_cycle(void)
{
	int nNumClut,nNumClut2;
	int alpha;

	for(nNumClut=1;nNumClut<=prot.DOF-1;++nNumClut)
	{
		for(nNumClut2=1;nNumClut2<=prot.DOF-1;++nNumClut2)
		{
			// 巨大行列の i,j 成分の計算
			Matrix[nNumClut][nNumClut2] = sub_calc_Matrix_cycle(nNumClut, nNumClut2);
		}
	}
}

// 巨大行列の i,j 成分の計算
double sub_calc_Matrix_cycle(int nNumClutI, int nNumClutJ)
{
	int alpha,i,j;
	int nNumClutMax;
	int nNumClut;

	double dummy = 0.0;

	nNumClutMax = nNumClutI;

	if (nNumClutMax < nNumClutJ)
	{
		nNumClutMax = nNumClutJ;
	}

	for (nNumClut=nNumClutMax;nNumClut<=prot.DOF-1;++nNumClut)
	{
		for (i=0;i<4;++i)
		{
			for (j=0;j<4;++j)
			{
				MatrixOfIJEle[i][j] = 0.0;
			}
		}

		// 巨大行列の i,j 成分のひとつの項の計算
		sub_sub_calc_Matrix_cycle(nNumClutI, nNumClutJ, nNumClut);
		for (i=0;i<4;++i)
		{
			dummy += MatrixOfIJEle[i][i];
		}
	}

	return dummy;
}

// 巨大行列の i,j 成分のひとつの項の計算
void sub_sub_calc_Matrix_cycle(int nNumClutI,
							   int nNumClutJ,
							   int nNumClutK){
  int i,j,k;
  
  double mat_dummy[4][4], matdotTone[4][4], matdotTtwo[4][4];

	for (i=0;i<4;++i)
	{
		for (j=0;j<4;++j)
		{
			matdotTone[i][j] = 0.0;
			matdotTtwo[i][j] = 0.0;
		}
	}

	calc_dot_Pseduo_TransMatrix(nNumClutJ, nNumClutK, matdotTone);
	calc_dot_Pseduo_TransMatrix(nNumClutI, nNumClutK, matdotTtwo);

	for (i=0;i<4;++i)
	{
		for (j=0;j<4;++j)
		{
			mat_dummy[i][j] = 0.0;
			MatrixOfIJEle[i][j] = 0.0;
		}
	}

	
	for (i=0;i<4;++i)
	{
		for (j=0;j<4;++j)
		{
			for (k=0;k<4;++k)
			{
				mat_dummy[i][j] += matdotTone[i][k]*clust[nNumClutK].PsedoInertia[k][j];
			}
		}
	}

	for (i=0;i<4;++i)
	{
		for (j=0;j<4;++j)
		{
			for (k=0;k<4;++k)
			{
				MatrixOfIJEle[i][j] += mat_dummy[i][k]*matdotTtwo[j][k];
			}
		}
	}

}

void calc_dot_Pseduo_TransMatrix(int nNumClutI,int nNumClutK,double mat[4][4]) {
  int i,j,k;
  int flag;
  int nNumClut,nNumClutdummy,nNumClutLast;
  double dummy_matrix[4][4], dummy_matrix2[4][4];

  flag=ON;
  if (clust[nNumClutI].join>0) {
    for (nNumClutdummy=nNumClutI;clust[nNumClutdummy].join!=clust[nNumClutI].join-1;++nNumClutdummy) {
      nNumClutLast = nNumClutdummy;
    }
    if (nNumClutLast < nNumClutK) {
      flag=OFF;
    }
  }
  for (i=0;i<4;++i) {
    for (j=0;j<4;++j) {
      mat[i][j] = 0.0;
    }
  }
  /**********************************************************************/
  /* printf("Mat:nNumClutI=%d,nNumClutK=%d 138\n",nNumClutI,nNumClutK); */
  /**********************************************************************/

  if (flag==ON) {
    for (i=0;i<4;++i){
      for (j=0;j<4;++j){
	dummy_matrix[i][j] = ident[i][j];
      }
    }
  
    for (nNumClut=1/*0*/;nNumClut<=nNumClutI;/* ++nNumClut */){
      if (clust[nNumClut].join <= clust[nNumClutI].join) {
	for (i=0;i<4;++i){
	  for (j=0;j<4;++j)	{
	    dummy_matrix2[i][j] = dummy_matrix[i][j];
	  }
	}
	/****************************************************/
        /* printf("Mat:nNumClut=%d 155\n",nNumClut);	    */
        /****************************************************/
  
	for (i=0;i<4;++i){
	  for (j=0;j<4;++j){
	    dummy_matrix[i][j] = 0.0;
	  }
	}
	/*******************************/
        /* printf("Mat:161\n");	       */
        /*******************************/

        /*************************************************************/
        /* printf("nNumClut=%d, nNumClutI=%d\n",nNumClut,nNumClutI); */
        /*************************************************************/
        
	for (i=0;i<4;++i){
	  for (j=0;j<4;++j){
	    for (k=0;k<4;++k){
	      dummy_matrix[i][j] += dummy_matrix2[i][k]*clust[nNumClut].PsedoTransMatrix[k][j];
	    }
	  }
	}
	/*******************************/
        /* printf("Mat:174\n");	       */
        /*******************************/
      }
      /************************/
      /* printf("Mat:177\n"); */
      /************************/
    
      if (clust[nNumClut].num_branch >1 ) {
	if (clust[nNumClut].nNumClutOfChild[1]-1 <= nNumClutI) {
	  nNumClut = clust[nNumClut].nNumClutOfChild[1]-1-1;
	}
      }
      ++nNumClut;
    }
    /************************/
    /* printf("Mat:185\n"); */
    /************************/
    
    for (i=0;i<4;++i) {
      for (j=0;j<4;++j) {
	dummy_matrix2[i][j] = dummy_matrix[i][j];
      }
    }
  
    for (i=0;i<4;++i) {
      for (j=0;j<4;++j) {
	dummy_matrix[i][j] = 0.0;
      }
    }
  
    for (i=0;i<4;++i) {
      for (j=0;j<4;++j) {
	for (k=0;k<4;++k) {
	  dummy_matrix[i][j] += dummy_matrix2[i][k]*delta_matrix[k][j];
	}
      }
    }
    /************************/
    /* printf("Mat:204\n"); */
    /************************/

    nNumClut=nNumClutI;
    if (clust[nNumClut].num_branch >1 ) {
      if (clust[nNumClut].nNumClutOfChild[1]-1 <= nNumClutK) {
	nNumClut = clust[nNumClut].nNumClutOfChild[1]-1-1;
      }
    }
    ++nNumClut;
    for (/*nNumClut=nNumClutI+1*/;nNumClut<=nNumClutK;/* ++nNumClut */) {
      if (clust[nNumClut].join <= clust[nNumClutK].join) {
	for (i=0;i<4;++i) {
	  for (j=0;j<4;++j)	{
	    dummy_matrix2[i][j] = dummy_matrix[i][j];
	  }
	}
	
	for (i=0;i<4;++i){
	  for (j=0;j<4;++j){
	    dummy_matrix[i][j] = 0.0;
	  }
	}
	/*******************************/
        /* printf("Mat:226\n");	       */
        /*******************************/
  
        /*************************************************************/
        /* printf("nNumClut=%d, nNumClutK=%d\n",nNumClut,nNumClutK); */
        /*************************************************************/
	
	for (i=0;i<4;++i){
	  for (j=0;j<4;++j)	{
	    for (k=0;k<4;++k) {
	      dummy_matrix[i][j] += dummy_matrix2[i][k]*clust[nNumClut].PsedoTransMatrix[k][j];
	    }
	  }
	}
	/*******************************/
        /* printf("Mat:239\n");	       */
        /*******************************/
      }
      if (clust[nNumClut].num_branch >1 ) {
	if (clust[nNumClut].nNumClutOfChild[1]-1 <= nNumClutK) {
	  nNumClut = clust[nNumClut].nNumClutOfChild[1]-1-1;
	}
      }
      /************************/
      /* printf("Mat:246\n"); */
      /************************/
      ++nNumClut;
    }
    
    for (i=0;i<4;++i) {
      for (j=0;j<4;++j){
	mat[i][j] = dummy_matrix[i][j];
      }
    }
    /************************/
    /* printf("Mat:256\n"); */
    /************************/
  }
}

