
double BAR_ite(double F,double *ene_U0_C0,double *ene_U0_C1, double *ene_U1_C0,double *ene_U1_C1, int n0, int n1, double criteria_BAR);
double MBAR_ite(double *F, double ***ene, int n_sim, int *n, double criteria_BAR, int MAXITE );
double MBAR_ite2(double **rou, double ***ene, int n_sim, int *n);

double MBAR_ite_high_speed_2(double *expF, double ***expene, int n_sim, int *n, double criteria_BAR, int MAXITE );

double *MBAR_AVE_oned(double *fene, double ***enek, int n_sim, int *n, double **od_data, double width,double *max,double *min, int *frame);

double *MBAR_AVE_oned_multi_temp(double *fene, double ***enek,double *betaene, int n_sim, int *n, double **od_data, double width,double *max,double *min, int *frame);

double *MBAR_AVE_oned_multi_temp_2(double *ef, double ***euk, int n_sim, int *n, double **od_data, double width,double *max,double *min, int *frame, int numk);

double MBAR_AVE_value_n_multi_temp(double *fene, double ***enek,double *betaene, int n_sim, int *n, double **data, double n_pow);

double *MBAR_ACM_histod(double *fene, double ***enek, int n_sim, int *n, double **od_data, double width,double *max,double *min, int *frame, double *W);
double *MBAR_ACM2_histod(double *fene, double ***enek, int n_sim, int *n, double **od_data, double width,double *max,double *min, int *frame, double *W, double *WT);

double *MBAR_ACM0(int n_sim, int *n, double *W, double *WT);
double *MBAR_ACM( int n_sim, int *n, double *W);
double *MBAR_ACM2(int n_sim, int *n, double *W, double *WT);

double *MBAR_setw(double *fene, double ***enek, int n_sim, int *n, double *W, double *WT);

double *MBAR_AVE_twod(double *fene, double ***enek, int n_sim, int *n, double ***td_data, double *width,double *max, double *min, int *frame);

double MBAR_AVEV_multi_temp(double *ef, double ***euk, int n_sim, int *n, double **od_data,int numk);

double MBAR_ite_high_speed(double *F, double ***ene, int n_sim, int *n, double criteria_BAR, int MAXITE );

double **MBAR_AVE_twod_multi_temp(double *ef, double ***euk, int n_sim, int *n, double ***td_data, double widthx, double widthy,double *maxx,double *minx,double *maxy,double *miny, int *framex, int *framey, int numk);

double MBAR_A(double *expF, double ***expene, int n_sim, int *n, double *A, double *dA);

double *MBAR_AVE_twod_wmaxmin(double *fene, double ***enek, int n_sim, int *n, double ***td_data, double *width,double *max, double *min, int *frame);
