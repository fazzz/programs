#ifndef INCLUDE_SBFF
#define INCLUDE_SBFF


struct potential_SBAA {
  // energy terms
  double p_t;

  double p_cnb_t,p_nnb_t;
  double p_d_t,p_a_t,p_b_t;

  double *p_cnb,*p_nnb;
  double *p_d,*p_a,*p_b;

  // parameters
  double *ALJ,*BLJ;
  double *cradii;
  double ALJnnb;
  double cradiinnb;

  double *BEQ,*AEQ,*DEQ;
  int numcnb,numnnb,*indexcnb,*indexnnb;

  double criteria;
  double EBB,ESC;
  double EI,EB,EA;
  double ecnb,ennb;
};

struct potential_SBAAMB {
  double p_t;

  struct potential_SBAA *ene1,*ene2;

};


int SBAAff_calcFFCNB(double *ALJ, 
		     double *BLJ,
		     double *p_cnb,
		     int numcnb, 
		     int *indexcnb,
		     int num_atom,
		     double *cord);

void SBAAff_set_CNB_PARM(double *cradii,
			 int numcnb,
			 double ecnb, 
			 double *ALJ, 
			 double *BLJ, 
			 int numatom);

int SBAAff_calcFFNNB(double ALJ,
		     double *p_nnb,
		     int numnnb, 
		     int *indexnnb,
		     int num_atom,
		     double *cord);

void SBAAff_set_NNB_PARM(double cradiinnb,
			 int numnnb,
			 double ennb,
			 double *ALJ
			 ,int numatom);

int SBAAff_calcDIHE(double *p_d,
		    double *cord,
		    double *DEQ,
		    double EBB,
		    double ESC,
		    double EI);

int SBAAff_calcANGLE(double *p_a,
		     double *cord, 
		     double EA, 
		     double *AEQ);

int SBAAff_calcBOND(double *p_b,
		    double *cord, 
		    double EB, 
		    double *BEQ);

int SBAAff_set_calcff(struct potential_SBAA *ene, 
		      double *cord,
		      int numatom);

double SBAAff_calcff(double *crd, 
		     int numatom,
		     struct potential_SBAA *ene);

double *SBAAff_make_nc_list(int *numnc,
			    double *cord, 
			    int numatom, 
			    double criteria);

int *SBAAff_make_nn_list(int *numnn, 
			 double *cord, 
			 int numatom, 
			 double criteria);

int SBAAff_get_eq_val(double *cord, 
		      double *BEQ, 
		      double *AEQ, 
		      double *DEQ);

int SBAAff_check_parameters(struct potential_SBAA *ene, 
			    int numatom);

int SBAAMBff_set_calcff(struct potential_SBAA *ene1,
			struct potential_SBAA *ene2,
			double *deltaV, 
			double *cord1,
			double *cord2,
			int numatom);

double SBAAMBff_calcff(double *crd,
		       int numatom,
		       struct potential_SBAA *ene1, 
		       struct potential_SBAA *ene2, 
		       double delta, 
		       double deltaV);

int SBAAff_set_protein_calcff(struct potential_SBAA *ene,
			      double *cord,
			      int numatom);

double *SBAAff_make_nc_protein_list( int *numnc, 
				     double *cord, 
				     int numatom, 
				     double criteria);

int SBAAMBff_set_protein_calcff(struct potential_SBAA *ene1,
				struct potential_SBAA *ene2,
				double *deltaV, 
				double *cord1,
				double *cord2,
				int numatom);

int *SBAAff_make_nn_protein_list(int *numnn, 
				 double *cord, 
				 int numatom, 
				 double criteria);

void SBAAff_set_parameters_protein_default(struct potential_SBAA *ene);

int within_3_neibor_res(int num1, int num2);

int resnum(int numatom);

#endif
