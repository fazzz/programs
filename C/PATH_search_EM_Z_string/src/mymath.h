
double Fermi(double x);
double calc_ave(int num, double *data);
double calc_msq(int num, double *data);
double calc_var(int num, double *data);
double calc_flu(int num, double *data);
double max_f(int num,double *data);
double min_f(int num,double *data);
void mac_min_f(int num,double *data, double *max, double *min);
double ts_normalize(double *timeseries, int nums, int numv, double *timeseriesnorm);
double *ts_normalize_od(double *timeseries,int numstep);
double ts_normalize_v2(double **timeseries, int nums, int numv, double **timeseriesnorm);
int ts_bl_ave(double **timeseries, int nts, int numv, int nbts, int nvts, double **timeseriesba);
double mm_len(double *vec);
int mm_outprod(double v1[3],double v2[3],double *v1x2);

double MYMATH_ave(double *data,int numstep);
double MYMATH_sqave(double *data,int numstep);
double MYMATH_var(double *data,int numstep);
