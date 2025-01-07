
#ifndef INCLUDE_GOLMAA
#define INCLUDE_GOLMAA

#define UNIT 4.184070*100.0

struct potential_GOLMAA {
  // energy terms
  double p_t;

  double p_natatt_t,p_repul_t;
  double p_d_t,p_a_t,p_b_t;
  double *p_natatt,*p_repul;
  double *p_d,*p_a,*p_b;

  // force terms
  double **f_t;

  double **f_natatt,**f_repul;
  double **f_d,**f_a,**f_b;

  // parameters
  double *ALJ_natatt,*BLJ_natatt;
  double ep_natatt;
  double *cradii_natatt;
  int  num_natatt,*index_natatt;
  int *ncmap;

  double ALJ_repul;
  double ep_repul;
  double cradii_repul;
  int num_repul,*index_repul;

  double FC_bond,FC_angle,FC_dihed;
  double *BEQ,*AEQ,*DEQ;
  int num_bond,*index_bond;
  int num_angl,*index_angl;
  int num_dihe,*index_dihe;

};

double GOLMAAff_calcff(double *crd, int numatom,struct potential_GOLMAA *ene,
		       int flagd, int flagnc, int flagnn,int **nb_matrix);

void GOLMAApote_calcDIHE(double *p_d,double *crd,double *DEQ,double FC_dihed);
//void GOLMAAforc_calcDIHE(double **f_d,double *crd,double *DEQ,double FC_dihed1,double FC_dihed2,int num_dihed,int *index_dihed,int numatom);

int GOLMAApote_calcNatAtt(double *p_natatt,double *crd,double *ALJ_natatt,double *BLJ_natatt,int num_natatt,int numatom,int *index_natatt);
int GOLMAApote_calcRepul(double *p_repul,double *crd,double ALJ_repul,int numatom);
int GOLMAAforc_calcNatAtt(double **f_natatt,double *crd,double *ALJ_natatt,double *BLJ_natatt,int num_natatt,int numatom,int *index_natatt);

int GOLMAAforc_calcRepul(double **f_repul,double *crd,double ALJ_repul,int numatom);

int GOLMAApoteforc_calcNatAttandRepul(double *p_natatt,double *p_repul,double **f_natatt,double **f_repul,double *crd,double *ALJ_natatt,double *BLJ_natatt,double e_natatt,double ALJ_repul,int num_natatt,int numatom,int *ncmap,int **nb_matrix);

double GOLMAA_calcTorque(double *Q,double *crd,double *DEQ,double FC_dihed,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a);

int *GOLMAAL_make_inpindex(int *inpnumA,int *indexclut,int numclut,int *nNumClutOfParent,int *terminal_atom_a,int *origin_atom_a);

void GOLMAA_out_formated(FILE *outputfile,struct potential_GOLMAA e,double KE,int i,double dt);

#endif
