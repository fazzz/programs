
#ifndef INCLUDE_MT_HS
#define INCLUDE_MT_HS

double ffLc_calcffandforce_HS(double *crd, int numatom,struct potential *ene,struct force *f,struct AmberParmL ap);

double ffLc_calcffandforce_nb_HS(double *ele, double *ALJ, double *BLJ, double *f,
				 int numnb, int *indexnb,
				 int num_atom,double *cord);

double ffLc_calcffandforce_14_HS(double *ele, double *ALJ, double *BLJ, double *f,
				 int numnb, int *indexnb,
				 int num_atom,double *cord);

double ffLc_calcDIHE_ffandforce_Cartesian_HS(int **PH, int **PA,
					     double *PN, double *PK, double *PHASE,
					     double *f_d,double *cord,
					     int numatom, int NPHIH, int MPHIA);

double ffLc_calcANGLE_ffandforce_Cartesian_HS(int **TH, int **TA, double *TK, double *TEQ,
					      double *f_a,double *cord,
					      int numatom, int NTHETH, int MTHETA);

double ffLc_calcBOND_ffandforce_Cartesian_HS(int **BH, int **BA,double *RK, double *REQ,
					     double *f_b,double *cord,
					     int numatom, int NBONH, int MBONA);

int ffLc_set_calcffandforce_HS(struct potential *ene, struct force *f,struct AmberParmL ap);

#endif

