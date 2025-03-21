#ifndef INCLUDE_FFLc
#define INCLUDE_FFLc

//#define UNIT 4.184070*100.0

#include "FFL.h"
#include "PTLb.h"

int ffLc_calcFFNB(double *ele, double *ALJ, double *BLJ,
		  double *p_e,double *p_LJ,
		  double *f_e,double *f_LJ,
		  int numnb, int *indexnb,
		  int num_atom,
		  double *cord,
		  int flagp, int flagf,struct AmberParmL ap);

int ffLc_calcFFNB_14(double *ele, double *ALJ, double *BLJ,
		    double *p_e,double *p_LJ,
		    double *f_e,double *f_LJ,
		    int numnb, int *indexnb,
		    int num_atom,
		    double *cord,
		    int flagp, int flagf,struct AmberParmL ap);

int ffLc_calcFFNB_wpep(double *ele, double *ALJ, double *BLJ,
		     double *p_e,double *p_LJ,
		     double *f_e,double *f_LJ,
		     int numnb, int *indexnb,
		     int num_atom,
		     double *cord,
		     double *u_1_5, double fact,
		     int flagp, int flagf,struct AmberParmL ap);

void ffLc_set_NB_PARM(double *ele, double *ALJ, double *BLJ, int numatom,struct AmberParmL ap);

int ffLc_set_NB_index(int *indexnb,int numnb, int numatom,struct AmberParmL ap);

int ffLc_set_numnb(struct AmberParmL ap);

int ffLc_calcDIHE(double *p_d,
		double *n_d,
		double *cord,
		int flagp, int flagf, int flaginp,struct AmberParmL ap);

int ffLc_calcDIHE_wpep(double *p_d,
		     double *n_d,
		     double *cord,
		     double *u_d,
		     int flagp, int flagf, int flaginp,struct AmberParmL ap);

int ffLc_calc_spe_type_DIHE(double *p,
			  double *n,
			  double angle,
			  int type,struct AmberParmL ap);

int ffLc_calcDIHE_force_Cartesian(double *f_d,double *cord,struct AmberParmL ap);

int ffLc_calcANGLE(double *p_a,double *cord,struct AmberParmL ap);

int ffLc_calcBOND(double *p_b,double *cord,struct AmberParmL ap);
int ffLc_calcBOND_force_Cartesian(double *f_b,double *cord,struct AmberParmL ap);

int ffLc_set_calcff(int numnb, int num14,FILE *inputfile, struct potential *ene,struct AmberParmL ap);

int ffLc_set_calcffsp(struct potential *ene,struct AmberParmL ap);

int ffLc_set_calcffandforce(struct potential *ene, struct force *f,struct AmberParmL ap);

double ffLc_calcff(double *crd, int numatom, struct potential *ene,struct AmberParmL ap);

double ffLc_calcffandforce(double *crd, int numatom, struct potential *ene,struct force *f,struct AmberParmL ap);

void ffLc_set_non_bonding_index_1(int *numindexnb, int *numindex14,struct AmberParmL ap);

void ffLc_set_non_bonding_index_2(int *gindexnb,int *gindex14,struct AmberParmL ap);

int ffLc_which_calc_nb(int atomi,int atomj,struct AmberParmL ap);

int ffLc_which_calc_14_nb(int atomi,int atomj,struct AmberParmL ap);

///////////////////////////////////////////////////////////////////

double ffLc_calcffandforce_woH(double *crd, int numatom,struct potential *ene,struct force *f,struct AmberParmL ap);

int ffLc_calcDIHE_woH(double *p_d,
		double *n_d,
		double *cord,
		    int flagp, int flagf, int flaginp,struct AmberParmL ap);

int ffLc_calcDIHE_force_Cartesian_woH(double *f_d,double *cord,struct AmberParmL ap);

int ffLc_calcANGLE_woH(double *p_a,double *cord,struct AmberParmL ap);

int ffLc_calcANGLE_force_Cartesian_woH(double *f_a,double *cord,struct AmberParmL ap);

int ffLc_calcBOND_woH(double *p_b,double *cord,struct AmberParmL ap);

int ffLc_calcBOND_force_Cartesian_woH(double *f_b,double *cord,struct AmberParmL ap);

int ffLc_calcFFNB_woH(double *ele, double *ALJ, double *BLJ,
		    double *p_e,double *p_LJ,
		    double *f_e,double *f_LJ,
		    int numnb, int *indexnb,
		    int num_atom,
		    double *cord,
		    int flagp, int flagf,struct AmberParmL ap);

int ffLc_calcFFNBvdW_woH( double *ALJ, double *BLJ,
			 double *p_LJ,
			 double *f_LJ,
			 int numnb, int *indexnb,
			 int num_atom,
			 double *cord,
			 int flagp, int flagf,struct AmberParmL ap);

int ffLc_calcBOND_woH_wFC100(double *p_b,double *cord,struct AmberParmL ap);

int ffLc_calcBOND_force_Cartesian_woH_eFC100(double *f_b,double *cord,struct AmberParmL ap);

///////////////////////////////////////////////////////////////////

int ffLc_calcFFNB_14_6(double *ele, double *ALJ, double *BLJ, // 0911
		     double *p_e,double *p_LJ,
		     double *f_e,double *f_LJ,
		     int numnb, int *indexnb,
		     int num_atom,
		     double *cord,
		     int flagp, int flagf,double scale,int numrepul,struct AmberParmL ap);

int ffLc_calcFFNB_wtune(double *ele, double *ALJ, double *BLJ,
		       double *p_e,double *p_LJ,
		       double *f_e,double *f_LJ,
		       int numnb, int *indexnb,
		       int num_atom,
		       double *cord,
		       int flagp, int flagf,
		       int *atom_tune_pairs, double *tune_val, int numtune,struct AmberParmL ap);

int ffLc_calcFFNB_wtuneb(double *ele, double *ALJ, double *BLJ,
			double *p_e,double *p_LJ,
			double *f_e,double *f_LJ,
			int numnb, int *indexnb,
			int num_atom,
			double *cord,
			int flagp, int flagf,
			int *atom_tune_pairs_es, int *atom_tune_pairs_LJ,
			double *tune_val_es, double *tune_val_LJ, 
			int numtune_es, int numtune_LJ,struct AmberParmL ap);

double ffLc_calcffandforce_w14tune(double *crd, int numatom,struct potential *ene,struct force *f, int *atom_tune_pairs, double *tune_val, int numtune,struct AmberParmL ap);

double ffLc_calcffandforce_w14tuneb(double *crd, int numatom,struct potential *ene,struct force *f, int *atom_tune_pairs_es, int *atom_tune_pairs_LJ, double *tune_val_es,double *tune_val_LJ, int numtune_es, int numtune_LJ,struct AmberParmL ap);

double ffLc_calcffandforce_14_6(double *crd, int numatom,struct potential *ene,struct force *f,
			      double scale,int numrepul,struct AmberParmL ap);

double ffLc_calcffandforce_14D_woH(double *crd, int numatom,struct potential *ene,struct force *f,struct AmberParmL ap);

double ffLc_calcffandforce_14DAB_woH(double *crd, int numatom,struct potential *ene,struct force *f,struct AmberParmL ap);

double ffLc_calcffandforce_14vdWDAB_woH(double *crd, int numatom,struct potential *ene,struct force *f,struct AmberParmL ap);

int ffLc_calcTorque_woH(double *Q,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a,struct AmberParmL ap);

int ffLc_calcTorque(double *Q,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a,struct AmberParmL ap);

double ffLc_calcTorque_wtune(double *Q,double *p_d,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a,int *atom_tune_pairs, double *tune_val_by_period, int numtune,struct AmberParmL ap);

double ffLc_calcTorque_wtuneb(double *Q,double *p_d,double *cord,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a,int *atom_tune_pairs, double *tune_val_V_n, double *tune_val_n_phase, int numtune, double pi,struct AmberParmL ap);

int *ffLc_make_inpindex(int *inpnumH,int *inpnumA,int *indexclut,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a,struct AmberParmL ap);

int *ffLc_make_inpindex_A(int *inpnumA,struct AmberParmL ap);

int ffLc_calcffandforce_speD(double *f_d,double *p_d,double *cord,int numspedihed,
			    int *atom1,int *atom2,int *atom3,int *atom4,struct AmberParmL ap);




void ffLc_out_formated(FILE *outputfile,struct potential e,double KE,double KEv,double PEv,double T,int i,double dt,struct AmberParmL ap);


#endif
