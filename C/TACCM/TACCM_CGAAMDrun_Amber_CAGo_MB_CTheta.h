
#ifndef INCLUDE_TA_MD_CE
#define INCLUDE_TA_MD_CE

double TACCM_CTheta_Amber_CAGo_MB(double *crd,int numatom,double *theta, 
				  int numdihe, int **pairs_dih_AA,
				  int numangl, int **pairs_ang_AA,
				  int numbond, int **pairs_bon_AA, 
				  double pi);

#endif
