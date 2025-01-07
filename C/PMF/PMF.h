

double *pmf_2dmap(double *data, int numdata,double width,double *maxx,double *maxy,double *minx,double *miny, int *framex,int *framey );

double *pmf_2dmap_rew(double *data,double *w, int numdata,double width,double *maxx,double *maxy,double *minx,double *miny,int *framex,int *framey,double temp);

double *pmf_1dmap(double *data, int numdata,double width,double *max,double *min,int *frame);

double *pmf_2dmap_wmaxmin(double *data, int numdata,double widthx,double widthy,double maxx,double maxy,double minx,double miny,int *framex,int *framey );

double *pmf_2dmap_wmaxmin_rew(double *data, double *w, int numdata,double widthx,double widthy,double maxx,double maxy,double minx,double miny,int *framex,int *framey, double T );
