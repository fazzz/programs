
double GearsConst5[6];
double GearsConst2[6];
double GearsConst4[5];
double Telar_Matrix[6][6];


void intg_set(void);
void intg_pc5(int flag,double prev[6],double corv[6], double dt,double acc);
void intg_pc4(int flag,double prev[5],double corv[5], double dt,double vel);
