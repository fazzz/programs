
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "ABAb.h"
#include "ABA_hosoku.h"
#include "preABAMD.h"

#include "EF.h"

double cp_clustdata(CLTh *clt1, CLTh *clt2){
  int i;

  (*clt2).origin_atom_a=(*clt1).origin_atom_a;
  (*clt2).terminal=(*clt1).terminal;
  (*clt2).num_atom_clust=(*clt1).num_atom_clust;
  (*clt2).num_branch=(*clt1).num_branch;
   
  for(i=0;i<(*clt2).num_branch;++i)
    (*clt2).terminal_atom_a[i]=0;
  //  (*clt2).terminal_atom_a=(int *)gcerealloc((*clt2).terminal_atom_a,sizeof(int)*(*clt1).num_branch);
  for(i=0;i<(*clt1).num_branch;++i)
    (*clt2).terminal_atom_a[i]=(*clt1).terminal_atom_a[i];
  (*clt2).nNumClutOfParent=(*clt1).nNumClutOfParent;
  
  for(i=0;i<(*clt2).num_branch;++i)
    (*clt2).nNumClutOfChild[i]=0;
  //  (*clt2).nNumClutOfChild=(int *)gcerealloc((*clt2).nNumClutOfChild,sizeof(int)*(*clt1).num_branch);
  for(i=0;i<(*clt1).num_branch;++i)
    (*clt2).nNumClutOfChild[i]=(*clt1).nNumClutOfChild[i];

  return 0.0;
}
