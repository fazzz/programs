#ifndef INCLUDE_ABA_M
#define INCLUDE_ABA_M

void ABAm_mainpass(ABI* abi,CLT *clt,double *Q,int numclut);

void ABAm_mainpass_TermOn(ABI* abi,CLT *clt,double *Q,int numclut,double *acc_Term,double *acc_Term2,double *vel_Term);

int ABAm_cPABI(double preABI[6][6],double corABI[6][6]);
void ABAm_cKG(double KG[6],double corABI[6][6]);
double ABAm_cD(double corABI[6][6]);
void ABAm_cCABI(double corABI[6][6],double preABI[4][6][6],double IM[6][6], double TM[4][6][6],int numbranch);

void ABAm_preBF(double preBF[6],double corBF[6],double KG[6], double eata);
double ABAm_eata(double T,double Corbf[6]);
double ABAm_nyu(double eata,double D);
void ABAm_corBF(double corBF[6],double preBF[4][6],double corABI[6][6],double Coracc[6],double Corfrc[6],double frc[6],double TM[4][6][6],int numbranch);

void ABAm_corBF_TERM(double corABI[6][6],double corBF[6],double *acc_Term,double *acc_Term2,double *vel_Term);

#endif
