
//#include "ParmTop.h"

int disulfid_count_numcys(void);

int disulfid_read_parm(int **atoms_b,double *K_bond,double *eq_bond,int *nbond,
		       int **atoms_a,double *K_angl,double *eq_angl,int *nangl,
		       int **atoms_d,double *V_dihe,double *n_dihe,double *the_dihe,int *ndihed
		       );

void disulfid_calc_pf(double *crd,double *p,double *f,
		      int **atoms_b,double *K_bond,double *eq_bond,int nbond,
		      int flagb,
		      int **atoms_a,double *K_angl,double *eq_angl,int nangl,
		      int flaga,
		      int **atoms_d,double *V_dihe,
		      double *n_dihe,double *the_dihe,int ndihed,
		      int flagd
		      );

void disulfid_calcBOND(double *crd,double *p_b,double *f_b,
		       int **atoms_b,double *K_bond,double *eq_bond,
		       int numbond);

void disulfid_calcANGL(double *crd,double *p_a,double *f_a,
		       int **atoms_a,double *K_angl,double *eq_angl,
		       int numangl);

void disulfid_calcDIHE(double *crd,double *p_d,double *f_d,
		       int **atoms_d,double *V_dihe,double *n_dihe,double *the_dihe,
		       int numdihed);

double disulfid_check_force_calc(double *crd,double *f,double dx,
				 int **atoms_b,double *K_bond,double *eq_bond,int nbond,
				 int flagb,
				 int **atoms_a,double *K_angl,double *eq_angl,int nangl,
				 int flaga,
				 int **atoms_d,double *V_dihe,
				 double *n_dihe,double *the_dihe,int ndihed,
				 int flagd);

void disulfid_calcBOND_check(double *crd,int numatom,double dx,
			     double *f_b,int **atoms_b,
			     double *K_bond,double *eq_bond,
			     int numbond);


void disulfid_calcANGL_check(double *crd,int numatom,double dx,
			     double *f_a,
			     int **atoms_a,double *K_angl,double *eq_angl,
			     int numangl);

void disulfid_calcDIHE_check(double *crd,int numatom,double dx,
			     double *f_d,
			     int **atoms_d,double *V_dihe,
			     double *n_dihe,double *the_dihe,
			     int numdihed);
