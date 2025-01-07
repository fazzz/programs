#ifndef INCLUDE_ABAb_sT
#define INCLUDE_ABAb_sT

void ABAbs_trans_mainpass(int stref,int endref, CLTb *clt,
			 double* PA1_t,double* PA2_t,double* PA21_t,double* bA1_t ,double* bA2_t,
			 double* PA1,double* PA2,double* PA21,double* bA1 ,double* bA2);
void ABAbs_trans_backpass();
void ABAbs_trans_b(double* b_trans,double *TMat,double *b);
void ABAbs_trans_P(double* P_trans,double* TMat,double* P);
void ABAbs_mak_transMat(double* transMat,CLTb *clt,int nstart, int nend);

#endif
