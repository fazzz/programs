//#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "const.h"
#include "dPCA.h"

void proj_dPCA(/*double *dihed_traj,*/double *dihed_ave,int time, int numdihed,
	       double covR[MAXNUMDIHED*MAXNUMDIHED], double w[MAXNUMDIHED])
{
  int i,j,jj;

  dihed_traj_trans=malloc(sizeof(double)*time*numdihed);
  for (i=0;i<time;++i)
    for (j=0;j<numdihed;++j)
      dihed_traj_trans[i*numdihed+j]=0.0;

  j=0;
  for (i=0;i<time;++i){
    for (j=0;j<2;++j){
	for (jj=0;jj<numdihed;++jj){
	  dihed_traj_trans[i*2+j]+=covR[jj*numdihed+j]*(dihed_traj[i*numdihed+jj]-dihed_ave[jj]);
	}
    }
  }

}
