
#ifndef INCLUDE_GOLM
#define INCLUDE_GOLM

#define UNIT 4.184070*100.0

struct potential_GOLM {
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

  double ALJ_repul;
  double ep_repul;
  double cradii_repul;
  int num_repul,*index_repul;

  double FC_bond,FC_angle,FC_dihed1,FC_dihed2;
  double *BEQ,*AEQ,*DEQ;
  int num_bond,*index_bond;
  int num_angl,*index_angl;
  int num_dihe,*index_dihe;

};

double GOLMff_calcff(double *crd, int numatom,struct potential_GOLM *ene,
		     int flagb,int flagc, int flagd, int flagnc, int flagnn);
int GOLMpote_calcBOND(double *p_b,double *cord, double FC_bond, double *BEQ,int num_bond,int *index_bond);
void GOLMforc_calcBOND(double **f_b,double *cord, double FC_bond, double *BEQ,int num_bond,int *index_bond, int numatom);
int GOLMpote_calcANGLE(double *p_a,double *crd,double FC_angle,double *AEQ,int num_angl,int *index_angl);
int GOLMforc_calcANGLE(double **f_a,double *crd,double FC_angle,double *AEQ,int num_angl,int *index_angl,int numatom);
void GOLMpote_calcDIHE(double *p_d,double *crd,double *DEQ,double FC_dihed1,double FC_dihed2,int num_dihed,int *index_dihed);
void GOLMforc_calcDIHE(double **f_d,double *crd,double *DEQ,double FC_dihed1,double FC_dihed2,int num_dihed,int *index_dihed,int numatom);

int GOLMpote_calcNatAtt(double *p_natatt,double *crd,double *ALJ_natatt,double *BLJ_natatt,double ep,int num_natatt,int numatom,int *index_natatt);
int GOLMpote_calcRepul(double *p_repul,double *crd,double ALJ_repul,int numatom);
int GOLMforc_calcNatAtt(double **f_natatt,double *crd,double *ALJ_natatt,double *BLJ_natatt,double ep,int num_natatt,int numatom,int *index_natatt);
int GOLMforc_calcRepul(double **f_repul,double *crd,double ALJ_repul,int numatom);

#endif
