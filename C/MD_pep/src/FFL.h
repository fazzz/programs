#ifndef INCLUDE_FF
#define INCLUDE_FF

#define UNIT 4.184070*100.0

struct ffL_parameters {
  double *ele,*ALJ,*BLJ;
  int numnb,num14,*indexnb,*index14;
  int numpara;
};

struct potential {
  double p_t;

  double p_e_t,p_LJ_t;
  double p_e_14_t,p_LJ_14_t;
  double p_d_t,p_a_t,p_b_t;

  double *p_e,*p_LJ;
  double *p_e_14,*p_LJ_14;
  double *p_d,*p_a,*p_b;

  struct ffL_parameters parm;
};

struct force {
  double *f_t;

  double *f_e,*f_LJ;
  double *f_e_14,*f_LJ_14;
  double *f_d,*f_a,*f_b;

};

int ffL_set_calcffandforce(struct potential *ene, struct force *f);

void ffL_set_NB_PARM(double *ele, double *ALJ, double *BLJ, int numatom);
void ffL_set_non_bonding_index_1(int *numindexnb, int *numindex14);
void ffL_set_non_bonding_index_2(int *gindexnb,int *gindex14);
int which_calc_nb(int atomi,int atomj);
int which_calc_14_nb(int atomi,int atomj);

double ffL_calcffandforce(double *crd, int numatom, struct potential *ene,struct force *f);

int ffL_calcFFNB(double *ele, double *ALJ, double *BLJ,
		double *p_e,double *p_LJ,
		double *f_e,double *f_LJ,
		int numnb, int *indexnb,
		int num_atom,
		double *cord,
		int flagp, int flagf);

int ffL_calcDIHE(double *p_d,
		double *n_d,
		double *cord,
		int flagp, int flagf, int flaginp);
int ffL_calcDIHE_force_Cartesian(double *f_d,double *cord);

int ffL_calcANGLE(double *p_a,double *cord);
int ffL_calcANGLE_force_Cartesian(double *f_a,double *cord);

double calcANGKE_force(double atomi[3],double atomj[3],double atomk[3],double kang,double ang_eq,double *f);

int ffL_calcBOND(double *p_b,double *cord);
int ffL_calcBOND_force_Cartesian(double *f_b,double *cord);

///////////////////////////////////////////////////////////////////

#endif
