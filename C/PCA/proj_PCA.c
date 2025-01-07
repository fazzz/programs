
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "PCA.h"

void proj_PCA(double *traj,double *traj_trans,int numstep,int numatom,
	      double *covR, double w[3*MAXNUMATOM])
{
  int i,j,jj,k,kk;
  double summass;
  double c;


  for (i=0;i<time;++i){
    //    for (j=0;j<numatom;++j){
    for (k=0;k<2/*3*/;++k){
	for (jj=0;jj<numatom;++jj){
	  for (kk=0;kk<3;++kk){
	    traj_trans[i*2/**numatom*3+j*3*/+k]+=covR[(jj*3+kk)*numatom*3+(j*3+k)]*traj[i*numatom*3+jj*3+kk]-crd_ave[jj*3+kk]);
	  }
	}
	  //      }
    }
  }

  //  test=fopen("pca_axis_all.txt","r");

  //  pca

}
