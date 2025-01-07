
#define ON  1
#define OFF 0
#define LOG 2

double *hist_mk2dhist(double *data, double widthx, double widthy , int numdata, double *maxx, double *maxy, double *minx, double *miny,int *framex, int *framey, int normflag);

double *hist_mk2dhist_ts(double *data, double widthx, double widthy , int numdata, double *maxx, double *maxy, double *minx, double *miny,int *framex, int *framey);

double *hist_mkhist(double *data, double width, int numdata, double *max, double *min, int *frame, int normflag);

double *hist_mkhist_wp(double *data, double width, int numdata, double *max, double *min, int *frame, int normflag, double *prob);

double *hist_mk2dhist_wmaxmin(double *data, double widthx, double widthy , int numdata, double maxx, double maxy, double minx, double miny,int *framex, int *framey, int normflag);
