
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "DCA.h"
#include "EF.h"

CLT *DCAp_clustscan(FILE *input,int *numclut){
  int i,j,x;
  CLT *clt;
  int *IndexOfABICycle;

  fscanf(input,"%d",numclut);

  clt=(CLT *)gcemalloc(sizeof(CLT)*(*numclut));
  for(i=0;i<(*numclut);++i) {
    clt[i].Coacc=(double *)gcemalloc(sizeof(double)*6);
    clt[i].Cofrc=(double *)gcemalloc(sizeof(double)*6);
  }

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

  return clt;
}
