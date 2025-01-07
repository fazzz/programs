
double ssl_covm(double *timeseries, int nums, int numv, double *COVM);
double ssl_ave(double *timeseries, int nums, int numv, double *AVE);
double ssl_normalize(double *timeseries, int nums, int numv, double *timeseriesnorm);
int ssl_subgraalg(double *S, double rou, double *W, int numv, double *beta);
double ssl_gralasso(double *S, double rou, int numv, double *W, double sigma, double lambda, double *l, double *omega);
int ssl_gralassomain(double *COVM, double rou, int numv, double *Lambda, double *Sigma );
double ssl_KLdiv(double *Lambda, double *Sigma, double *Lambdap, double *Sigmap, int numv,  double *KLdiv, int topten[10], double vtopten[10]);
