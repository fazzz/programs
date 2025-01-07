#ifndef INCLUDE_FF
#define INCLUDE_FF

#define UNIT 4.184070*100.0

struct ff_parameters {
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

  struct ff_parameters parm;
};

struct force {
  double *f_t;

  double *f_e,*f_LJ;
  double *f_e_14,*f_LJ_14;
  double *f_d,*f_a,*f_b;

};


int ff_calcFFNB(double *ele, double *ALJ, double *BLJ,
		double *p_e,double *p_LJ,
		double *f_e,double *f_LJ,
		int numnb, int *indexnb,
		int num_atom,
		double *cord,
		int flagp, int flagf);

int ff_calcFFNB_wpep(double *ele, double *ALJ, double *BLJ,
		     double *p_e,double *p_LJ,
		     double *f_e,double *f_LJ,
		     int numnb, int *indexnb,
		     int num_atom,
		     double *cord,
		     double *u_1_5, double fact,
		     int flagp, int flagf);

void ff_set_NB_PARM(double *ele, double *ALJ, double *BLJ, int numatom);

int ff_set_NB_index(int *indexnb,int numnb, int numatom);

int ff_set_numnb(void);

int ff_calcDIHE(double *p_d,
		double *n_d,
		double *cord,
		int flagp, int flagf, int flaginp);

int ff_calcDIHE_wpep(double *p_d,
		     double *n_d,
		     double *cord,
		     double *u_d,
		     int flagp, int flagf, int flaginp);

int ff_calc_spe_type_DIHE(double *p,
			  double *n,
			  double angle,
			  int type);

int ff_calcDIHE_force_Cartesian(double *f_d,double *cord);

int ff_calcANGLE(double *p_a,double *cord);

int ff_calcBOND(double *p_b,double *cord);
int ff_calcBOND_force_Cartesian(double *f_b,double *cord);

int ff_set_calcff(int numnb, int num14,FILE *inputfile, struct potential *ene);

int ff_set_calcffsp(struct potential *ene);

int ff_set_calcffandforce(struct potential *ene, struct force *f);

double ff_calcff(double *crd, int numatom, struct potential *ene);

double ff_calcffandforce(double *crd, int numatom, struct potential *ene,struct force *f);

void ff_set_non_bonding_index_1(int *numindexnb, int *numindex14);

void ff_set_non_bonding_index_2(int *gindexnb,int *gindex14);

//int which_calc_nb(int atomi,int atomj); // 2014-09-19

//int which_calc_14_nb(int atomi,int atomj);

///////////////////////////////////////////////////////////////////

double ff_calcffandforce_woH(double *crd, int numatom,struct potential *ene,struct force *f);

int ff_calcDIHE_woH(double *p_d,
		double *n_d,
		double *cord,
		    int flagp, int flagf, int flaginp);

int ff_calcDIHE_force_Cartesian_woH(double *f_d,double *cord);

int ff_calcANGLE_woH(double *p_a,double *cord);

int ff_calcANGLE_force_Cartesian_woH(double *f_a,double *cord);

int ff_calcBOND_woH(double *p_b,double *cord);

int ff_calcBOND_force_Cartesian_woH(double *f_b,double *cord);

int ff_calcFFNB_woH(double *ele, double *ALJ, double *BLJ,
		    double *p_e,double *p_LJ,
		    double *f_e,double *f_LJ,
		    int numnb, int *indexnb,
		    int num_atom,
		    double *cord,
		    int flagp, int flagf);

int ff_calcBOND_woH_wFC100(double *p_b,double *cord);

int ff_calcBOND_force_Cartesian_woH_eFC100(double *f_b,double *cord);

///////////////////////////////////////////////////////////////////

int ff_calcFFNB_14_6(double *ele, double *ALJ, double *BLJ, // 0911
		     double *p_e,double *p_LJ,
		     double *f_e,double *f_LJ,
		     int numnb, int *indexnb,
		     int num_atom,
		     double *cord,
		     int flagp, int flagf,double scale,int numrepul);

double ff_calcffandforce_14_6(double *crd, int numatom,struct potential *ene,struct force *f,
			      double scale,int numrepul);

int ff_calcTorque(double *Q,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a);
int *make_inpindex(int *inpnumH,int *inpnumA,int *indexclut,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a);

#endif
