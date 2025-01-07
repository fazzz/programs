#ifndef INCLUDE_ABA_sT
#define INCLUDE_ABA_sT

void ABAs_trans_mainpass(int stref,int endref, CLT *clt,
			 double* PA1_t,double* PA2_t,double* PA21_t,double* bA1_t ,double* bA2_t,
			 double* PA1,double* PA2,double* PA21,double* bA1 ,double* bA2);
void ABAs_trans_backpass();
void ABAs_trans_b(double* b_trans,double *TMat,double *b);
void ABAs_trans_P(double* P_trans,double* TMat,double* P);
void ABAs_mak_transMat(double* transMat,CLT *clt,int nstart, int nend);

#endif
