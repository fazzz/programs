
#ifndef INCLUDE_GOLMAA_DB
#define INCLUDE_GOLMAA_DB

struct potential_GOLMAA_dbasin {
  struct potential_GOLMAA L1;
  struct potential_GOLMAA L2;

  int **nb_matrix1,**nb_matrix2;

  double p_t;
  double **f_t;
};

double GOLMAAff_dbasin_calcff_ratio(double *crd, int numatom,struct potential_GOLMAA_dbasin *ene,double delta, double deltaV);
double GOLMAA_dbasin_calcTorque(double *Q,double *crd,struct potential_GOLMAA_dbasin *ene,int numclut,int *nNumClutOfParent,int *terminal,int *origin,double delta, double deltaV);

int GOLMAAff_dbasin_set_calcff(struct potential_GOLMAA_dbasin *ene, double *refcrd1, double *refcrd2,int numatom,double R_C_D);

#endif
