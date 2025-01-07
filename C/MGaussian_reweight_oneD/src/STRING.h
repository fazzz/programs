
double z_string_sample(double *p_stat/*,double *p_evoluted,double *f*/,int numpoint,double dt);
// double z_string_sample(double *p_stat,double *p_evoluted,double *f,int numpoint,double dt);
// double z_string_sample(double *p,double *p_evoluted,double *force,int numpoint,double dt);
// double z_string_sample(double *p,int numpoint,double dt);
double z_string_evolution(double *p, double *p_evoluted, double *fe, int numpoint,double dt, int numdimension);
double z_string_interpolation(double *p, double *p_evoluted,int numpoint, int numdimension);
double z_string_Cartesian(double *path,double *path_evoluted,double *fe,int numpoint,int numatom,double dt);
// double z_string_Cartesian(double *path,int numpoint,int numatom,double dt);

// double z_string_FASYS_calcpote(double **q);
double z_string_FASYS_calcpote(double **q, double kd[2],double n[2]);
// double z_string_FASYS(double *path,int numpoint,double dt, double *v);
double z_string_FASYS(double *path,int numpoint,double dt, double *v, double kd[2], double n[2]);
// double z_string_FASYS(double *path,int numpoint,double dt);
// double z_string_FASYS_calcforce(double **q,double **f);
//double z_string_FASYS_calcforce(double **q,double **f, double kd[2], double n[2]);
double z_string_FASYS_calcforce(double **q,double **f, double kd[2], double n[2], double **fb, double **fa, double **fd);
