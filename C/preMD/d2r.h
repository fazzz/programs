
#ifndef INCLUDE_d2r
#define INCLUDE_d2r

#define ON 1
#define OFF 0

#define ATOM 0
#define RESIDUE 1

#define start 0
#define fin   1

int *readd2rinput(FILE *inputfile, int *numres, int specifymove);

void writed2routput(FILE *clustfileout,CLTh *clt, int numclut);

int **res2atom_BB(int *resid, int numatom, int numres, int *numd2r);

int **res2atom_SC(int *resid, int numres, int *numd2r,CLTh *clt, int numclut);

double trans_d2r(CLTh *clt, int *numclut, int **d2r, int numd2r);

#endif
