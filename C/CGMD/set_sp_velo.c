#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "physics.h"
#include "MD.h"

// spatial velocity �̐ݒ�̕⏕���s���֐�
void sub_sub_set_sp_velo(int nNumClut,
	                     int nNumClutminusone);

// spatial velocity �̐ݒ���s���֐�
void set_sp_velo(int nNumClut, int nNumClutOrigBranch){
  int i;
  int nNumClutOfParent;

  nNumClutOfParent = clust[nNumClut].nNumClutOfParent-1;

  // �N���X�^�����[�̂Ƃ�
  if (nNumClut == 0)	{
    if (TermMoveMode==ON) {
      for(i=0;i<6;++i) {
	clust[0].sp_velo[i] = vel_Term[i];
      }
    }
    else {
      for(i=0;i<6;++i) {
	clust[0].sp_velo[i] = 0.0;
      }
    }
  }
  else {
    sub_set_sp_velo(nNumClut, nNumClutOfParent);
  }
}

// spatial velocity �̐ݒ�̕⏕���s���֐�_1
void sub_set_sp_velo(int nNumClut, int nNumClutminusone) {
  int i,j,k;
  int alpha;
  int alpha2;
  int nNumFreedom;

  FILE *test;
  // spatial velocity �̏��������s��
  for(i=0;i<6;++i) {
    clust[nNumClut].sp_velo[i] = 0.0;
  }

  // spatial velocity �̐ݒ���s��_1
  for(i=0;i<6;++i) {
    for(j=0;j<6;++j) {
      // �ϊ��s��̓]�u�s����悶��
      clust[nNumClut].sp_velo[i]
	+=  clust[nNumClut].TransMatrix[0][j][i]*clust[nNumClutminusone].sp_velo[j];
    }
  }

  // ���R�x�̑I��
  nNumFreedom = ABI[nNumClut].hingmat-1;

  // spatial velocity �̐ݒ���s��_2
  clust[nNumClut].sp_velo[nNumFreedom]+= clust[nNumClut].ddihedang[0];

//	test = fopen("sp_velo.out", "a");
//	fprintf(test," %e \n",  clust[nNumClut].ddihedang[0]/**1.0e20*/);
//	fprintf(test,"sp_velo[%d] = ", nNumClut);
//	for (i=0;i<6;++i)
//	{
//		fprintf(test," %e \n",  clust[nNumClut].sp_velo[i]/**1.0e20*/);
//	}
//	fprintf(test, "\n");
//	fclose(test);

}

