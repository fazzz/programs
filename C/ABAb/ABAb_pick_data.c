
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABAb.h"

CLTb *ABAbp_clustscan(FILE *input,int *numclut){
  int i,j,x;
  CLTb *clt;
  int *IndexOfABIbCycle;

  fscanf(input,"%d",numclut);

  clt=(CLTb *)gcemalloc(sizeof(CLTb)*(*numclut));

  for(i=0;i<(*numclut);++i) fscanf(input,"%d",&clt[i].origin_atom_a);
  for(i=0;i<(*numclut);++i) fscanf(input,"%d",&clt[i].terminal);
  for(i=0;i<(*numclut);++i) fscanf(input,"%d",&clt[i].num_atom_clust);
  for(i=0;i<(*numclut);++i) {
    fscanf(input,"%d",&clt[i].num_branch);
    clt[i].terminal_atom_a=(int *)gcemalloc(sizeof(int)*clt[i].num_branch);
    clt[i].nNumClutOfChild=(int *)gcemalloc(sizeof(int)*clt[i].num_branch);
  }
  for(i=0;i<(*numclut);++i) fscanf(input,"%d",&x);
  for(i=0;i<(*numclut);++i) 
    for(j=0;j<clt[i].num_branch;++j)
      fscanf(input,"%d",&(clt[i].terminal_atom_a[j]));
  for (i=0;i<(*numclut);++i) fscanf(input, "%d", &clt[i].nNumClutOfParent);
  for (i=0;i<(*numclut);++i)
    for(j=0;j<clt[i].num_branch;++j)
      fscanf(input, "%d", &clt[i].nNumClutOfChild[j]);
  IndexOfABIbCycle=(int *)gcemalloc(sizeof(int)*(*numclut));
  for (i=0;i<(*numclut);++i) fscanf(input, "%d", &IndexOfABIbCycle[i]);

  for(i=0;i<(*numclut);++i) {
    clt[i].frc_e=(double *)gcemalloc(sizeof(double)*clt[i].num_atom_clust*3);
    clt[i].frc_LJ=(double *)gcemalloc(sizeof(double)*clt[i].num_atom_clust*3);
    clt[i].frc_1_4_e=(double *)gcemalloc(sizeof(double)*clt[i].num_atom_clust*3);
    clt[i].frc_1_4_LJ=(double *)gcemalloc(sizeof(double)*clt[i].num_atom_clust*3);
  }

  return clt;
}
