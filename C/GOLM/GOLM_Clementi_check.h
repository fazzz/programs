
#ifndef INCLUDE_GOLM_check_set
#define INCLUDE_GOLM_check_set

double GOLM_Clementi_calcff_check(char *inputfilename,char *reffilename,char *parmfilename,
				  int numspatom,double dx,
				  double f_natatt[3],double f_repul[3],double f_d[3],double f_a[3],double f_b[3],
				  double f_d1[4][3],double f_d2[4][3], int nums,
				  double f_natatt1[2][3],double f_natatt2[4][3], int numa);


double GOLM_Clementi_ff_calcDIHE_s(double *crd,int numatom,double f_d[4][3],double Kd1,double Kd2,double *dih_equ, int nums);

double GOLM_Clementi_ff_calcff_natatt_s(double *crd, int numatom,int *index_numatt,int numatt,
					double *ALJ_natatt,double *BLJ_natatt,double ep_natatt,
					double p_natatt,double f_natatt[2][3],int nums);

#endif

