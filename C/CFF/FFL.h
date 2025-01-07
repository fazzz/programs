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


int ffL_calcFFNB(double *ele, double *ALJ, double *BLJ,
		double *p_e,double *p_LJ,
		double *f_e,double *f_LJ,
		int numnb, int *indexnb,
		int num_atom,
		double *cord,
		int flagp, int flagf);

int ffL_calcFFNB_14(double *ele, double *ALJ, double *BLJ,
		    double *p_e,double *p_LJ,
		    double *f_e,double *f_LJ,
		    int numnb, int *indexnb,
		    int num_atom,
		    double *cord,
		    int flagp, int flagf);

int ffL_calcFFNB_wpep(double *ele, double *ALJ, double *BLJ,
		     double *p_e,double *p_LJ,
		     double *f_e,double *f_LJ,
		     int numnb, int *indexnb,
		     int num_atom,
		     double *cord,
		     double *u_1_5, double fact,
		     int flagp, int flagf);

void ffL_set_NB_PARM(double *ele, double *ALJ, double *BLJ, int numatom);

int ffL_set_NB_index(int *indexnb,int numnb, int numatom);

int ffL_set_numnb(void);

int ffL_calcDIHE(double *p_d,
		double *n_d,
		double *cord,
		int flagp, int flagf, int flaginp);

int ffL_calcDIHE_wpep(double *p_d,
		     double *n_d,
		     double *cord,
		     double *u_d,
		     int flagp, int flagf, int flaginp);

int ffL_calc_spe_type_DIHE(double *p,
			  double *n,
			  double angle,
			  int type);

int ffL_calcDIHE_force_Cartesian(double *f_d,double *cord);

int ffL_calcANGLE(double *p_a,double *cord);

int ffL_calcBOND(double *p_b,double *cord);
int ffL_calcBOND_force_Cartesian(double *f_b,double *cord);

int ffL_set_calcff(int numnb, int num14,FILE *inputfile, struct potential *ene);

int ffL_set_calcffsp(struct potential *ene);

int ffL_set_calcffandforce(struct potential *ene, struct force *f);

double ffL_calcff(double *crd, int numatom, struct potential *ene);

double ffL_calcffandforce(double *crd, int numatom, struct potential *ene,struct force *f);

void ffL_set_non_bonding_index_1(int *numindexnb, int *numindex14);

void ffL_set_non_bonding_index_2(int *gindexnb,int *gindex14);

int which_calc_nb(int atomi,int atomj);

int which_calc_14_nb(int atomi,int atomj);

///////////////////////////////////////////////////////////////////

double ffL_calcffandforce_woH(double *crd, int numatom,struct potential *ene,struct force *f);

int ffL_calcDIHE_woH(double *p_d,
		double *n_d,
		double *cord,
		    int flagp, int flagf, int flaginp);

int ffL_calcDIHE_force_Cartesian_woH(double *f_d,double *cord);

int ffL_calcANGLE_woH(double *p_a,double *cord);

int ffL_calcANGLE_force_Cartesian_woH(double *f_a,double *cord);

int ffL_calcBOND_woH(double *p_b,double *cord);

int ffL_calcBOND_force_Cartesian_woH(double *f_b,double *cord);

int ffL_calcFFNB_woH(double *ele, double *ALJ, double *BLJ,
		    double *p_e,double *p_LJ,
		    double *f_e,double *f_LJ,
		    int numnb, int *indexnb,
		    int num_atom,
		    double *cord,
		    int flagp, int flagf);

int ffL_calcFFNBvdW_woH( double *ALJ, double *BLJ,
			 double *p_LJ,
			 double *f_LJ,
			 int numnb, int *indexnb,
			 int num_atom,
			 double *cord,
			 int flagp, int flagf);

int ffL_calcBOND_woH_wFC100(double *p_b,double *cord);

int ffL_calcBOND_force_Cartesian_woH_eFC100(double *f_b,double *cord);

///////////////////////////////////////////////////////////////////

int ffL_calcFFNB_14_6(double *ele, double *ALJ, double *BLJ, // 0911
		     double *p_e,double *p_LJ,
		     double *f_e,double *f_LJ,
		     int numnb, int *indexnb,
		     int num_atom,
		     double *cord,
		     int flagp, int flagf,double scale,int numrepul);

int ffL_calcFFNB_wtune(double *ele, double *ALJ, double *BLJ,
		       double *p_e,double *p_LJ,
		       double *f_e,double *f_LJ,
		       int numnb, int *indexnb,
		       int num_atom,
		       double *cord,
		       int flagp, int flagf,
		       int *atom_tune_pairs, double *tune_val, int numtune);

int ffL_calcFFNB_wtuneb(double *ele, double *ALJ, double *BLJ,
			double *p_e,double *p_LJ,
			double *f_e,double *f_LJ,
			int numnb, int *indexnb,
			int num_atom,
			double *cord,
			int flagp, int flagf,
			int *atom_tune_pairs_es, int *atom_tune_pairs_LJ,
			double *tune_val_es, double *tune_val_LJ, 
			int numtune_es, int numtune_LJ);

double ffL_calcffandforce_w14tune(double *crd, int numatom,struct potential *ene,struct force *f, int *atom_tune_pairs, double *tune_val, int numtune);

double ffL_calcffandforce_w14tuneb(double *crd, int numatom,struct potential *ene,struct force *f, int *atom_tune_pairs_es, int *atom_tune_pairs_LJ, double *tune_val_es,double *tune_val_LJ, int numtune_es, int numtune_LJ);

double ffL_calcffandforce_14_6(double *crd, int numatom,struct potential *ene,struct force *f,
			      double scale,int numrepul);

double ffL_calcffandforce_14D_woH(double *crd, int numatom,struct potential *ene,struct force *f);

double ffL_calcffandforce_14DAB_woH(double *crd, int numatom,struct potential *ene,struct force *f);

double ffL_calcffandforce_14vdWDAB_woH(double *crd, int numatom,struct potential *ene,struct force *f);

int ffL_calcTorque_woH(double *Q,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a);

int ffL_calcTorque(double *Q,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a);

double ffL_calcTorque_wtune(double *Q,double *p_d,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a,int *atom_tune_pairs, double *tune_val_by_period, int numtune);

double ffL_calcTorque_wtuneb(double *Q,double *p_d,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a,int *atom_tune_pairs, double *tune_val_V_n, double *tune_val_n_phase, int numtune, double pi);

int *ffL_make_inpindex(int *inpnumH,int *inpnumA,int *indexclut,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a);

int *ffL_make_inpindex_A(int *inpnumA);

int ffL_calcffandforce_speD(double *f_d,double *p_d,double *cord,int numspedihed,
			    int *atom1,int *atom2,int *atom3,int *atom4);




void ffL_out_formated(FILE *outputfile,struct potential e,double KE,double KEv,double PEv,double T,int i,double dt);


#endif
