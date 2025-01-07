
#ifndef INCLUDE_d2rc
#define INCLUDE_d2rc

#define ON 1
#define OFF 0

#define ATOM 0
#define RESIDUE 1

#define start 0
#define fin   1

int *readd2rinput(FILE *inputfile, int *numres, int specifymove);

void writed2routput(FILE *clustfileout,CLTb *clt, int numclut);

int **res2atom_BBc(int *resid, int numatom, int numres, int *numd2r);

int **res2atom_SCc(int *resid, int numres, int *numd2r,CLTb *clt, int numclut);

double trans_d2rc(CLTb *clt, int *numclut, int **d2r, int numd2r);

#endif
