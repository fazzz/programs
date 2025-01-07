#ifndef INCLUDE_ABAb_M
#define INCLUDE_ABAb_M

void ABAbm_mainpass(ABIb* abi,CLTb *clt,double *Q,int numclut);

void ABAbm_mainpass_TermOn(ABIb* abi,CLTb *clt,double *Q,int numclut,double *acc_Term,double *acc_Term2,double *vel_Term);

int ABAbm_cPABIb(double preABIb[6][6],double corABIb[6][6]);
void ABAbm_cKG(double KG[6],double corABIb[6][6]);
double ABAbm_cD(double corABIb[6][6]);
void ABAbm_cCABIb(double corABIb[6][6],double ***preABIb/*[10][6][6]*/,double IM[6][6], double ***TM/*[10][6][6]*/,int numbranch);

void ABAbm_preBF(double preBF[6],double corBF[6],double KG[6], double eata);
double ABAbm_eata(double T,double Corbf[6]);
double ABAbm_nyu(double eata,double D);
void ABAbm_corBF(double corBF[6],double **preBF/*[10][6]*/,double corABIb[6][6],double Coracc[6],double Corfrc[6],double frc[6],double ***TM/*[10][6][6]*/,int numbranch);

void ABAbm_corBF_TERM(double corABIb[6][6],double corBF[6],double *acc_Term,double *acc_Term2,double *vel_Term);

#endif
