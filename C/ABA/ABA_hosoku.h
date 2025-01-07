#ifndef INCLUDE_ABA_h
#define INCLUDE_ABA_h

struct clustdata_h{
  int origin_atom_a;
  int terminal_atom_a[/*4*/100]; // 2012-01-27
  int terminal;
  int num_clust;
  int num_atom_clust;
  int num_branch;
  int join;

  int nNumClutOfParent;
  int nNumClutOfChild[/*4*/100]; // 2012-01-27
  
};

typedef struct clustdata_h CLTh;

#define TERM 0

int ABAh_setJoin(CLTh *clt,int nNumClut);
CLTh *ABAhp_clustscan(FILE *input,int *numclut);

#endif

