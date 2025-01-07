
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ABA_hosoku.h"

int ABAh_setJoin(CLTh *clt,int nNumClut) {
  int i,j,k;
  int nNumParent;

  nNumParent = clt[nNumClut].nNumClutOfParent-1;

  if ( clt[nNumParent].num_branch > 1 && nNumParent > -1) {
    for (i=0;i<clt[nNumParent].num_branch;++i)
      if (nNumClut==clt[nNumParent].nNumClutOfChild[i]-1)
	break;

    if (i!=clt[nNumParent].num_branch-1)
      clt[nNumClut].join = clt[nNumParent].join+1/*+i*/;
    else
      clt[nNumClut].join = clt[nNumParent].join;
  }
  else if (nNumParent > -1)
    clt[nNumClut].join = clt[nNumParent].join;

  /***************************************/
  /* clt[nNumClut].join = joinflag;	 */
  /* 					 */
  /* if (clt[nNumClut].num_branch > 1) { */
  /*   joinflag +=1;			 */
  /* }					 */
  /* if (clt[nNumClut].terminal == 0) {	 */
  /*   joinflag -=1;			 */
  /* }					 */
  /***************************************/

  //  return joinflag;
}

CLTh *ABAhp_clustscan(FILE *input,int *numclut){
  int i,j,x;
  CLTh *clt;
  int *IndexOfABIbCycle;

  fscanf(input,"%d",numclut);

  clt=(CLTh *)gcemalloc(sizeof(CLTh)*(*numclut));

  for(i=0;i<(*numclut);++i) fscanf(input,"%d",&clt[i].origin_atom_a);
  for(i=0;i<(*numclut);++i) fscanf(input,"%d",&clt[i].terminal);
  for(i=0;i<(*numclut);++i) fscanf(input,"%d",&clt[i].num_atom_clust);
  for(i=0;i<(*numclut);++i) {
    fscanf(input,"%d",&clt[i].num_branch);
    //    clt[i].terminal_atom_a=(int *)gcemalloc(sizeof(int)*clt[i].num_branch);
    //    clt[i].nNumClutOfChild=(int *)gcemalloc(sizeof(int)*clt[i].num_branch);
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

  return clt;
}
