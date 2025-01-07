
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "ABA.h" // 2014-06-18
#include "ABAb.h"  // 2014-06-18
#include "EF.h" // 2014-06-17

CLT *ABAp_clustscan(FILE *input,int *numclut){ // 2014-06-20
//int ABAp_clustscan(FILE *input,int *numclut,CLT * clt){ // 2014-06-20
  int i,j,x;
  CLT *clt; // 2014-06-20
  int *IndexOfABICycle;

  fscanf(input,"%d",numclut);

  clt=(CLT *)gcemalloc(sizeof(CLT)*(*numclut));

  for(i=0;i<(*numclut);++i) fscanf(input,"%d",&clt[i].origin_atom_a);
  for(i=0;i<(*numclut);++i) fscanf(input,"%d",&clt[i].terminal);
  for(i=0;i<(*numclut);++i) fscanf(input,"%d",&clt[i].num_atom_clust);
  for(i=0;i<(*numclut);++i) fscanf(input,"%d",&clt[i].num_branch);
  for(i=0;i<(*numclut);++i) fscanf(input,"%d",&x);
  for(i=0;i<(*numclut);++i) 
    for(j=0;j<clt[i].num_branch;++j)
      fscanf(input,"%d",&clt[i].terminal_atom_a[j]);
  for (i=0;i<(*numclut);++i) fscanf(input, "%d", &clt[i].nNumClutOfParent);
  for (i=0;i<(*numclut);++i)
    for(j=0;j<clt[i].num_branch;++j)
      fscanf(input, "%d", &clt[i].nNumClutOfChild[j]);
  IndexOfABICycle=(int *)gcemalloc(sizeof(int)*(*numclut));
  for (i=0;i<(*numclut);++i) fscanf(input, "%d", &IndexOfABICycle[i]);

  for(i=0;i<(*numclut);++i) {
    clt[i].frc_e=(double *)gcemalloc(sizeof(double)*clt[i].num_atom_clust*3);
    clt[i].frc_LJ=(double *)gcemalloc(sizeof(double)*clt[i].num_atom_clust*3);
    clt[i].frc_1_4_e=(double *)gcemalloc(sizeof(double)*clt[i].num_atom_clust*3);
    clt[i].frc_1_4_LJ=(double *)gcemalloc(sizeof(double)*clt[i].num_atom_clust*3);
  }

  return clt; // 2014-06-20
  //  return 0;       // 2014-06-20
}
