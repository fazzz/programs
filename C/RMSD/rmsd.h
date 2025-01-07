#define AA 0
#define CA 1
#define HV 2

double rmsd_qcp(double *coordA,double *coordB,int numatom,int MODE);
void mmult(double *coordA, double *coordB, double matrix[3][3], int numatom, double *mass);
//void mmult(double *coordA, double *coordB, double matrix[3][3], int numatom);
void fomKmat(double Kmat[4][4], double mat[3][3]);
void fomPolynominal(/*double C0, double C1, double C2,*/ double mat[3][3], double Kmat[4][4]);
double Newton_Rapson(/*double C0, double C1, double C2,*/double lambda_ini, double det);
double CalcG(double *coordA,int numatom);
void transCentroid(double *coordA, double *coordB, int numatom, double *mass);

double InnerProduct(double *A, double **coords1, double **coords2, const int len, const double *weight);
int FastCalcRMSDAndRotation(double *rot, double *A, double *rmsd, double E0, int len, double minScore);
void CenterCoords(double **coords, const int len, const double *weight);
double CalcRMSDRotationalMatrix(double **coords1, double **coords2, const int len, double *rot, const double *weight);

