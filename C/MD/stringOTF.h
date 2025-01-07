#ifndef INCLUDE_STOF
#define INCLUDE_STOF

#include "FFL.h"

double String_Propagetor_Iso_FASYS(double *crd,double *vel,double *mass,int numatom,double IsoCoff,double dt,double *KE,double *PE,struct potential e,struct force f,double *z,double kappa,int numcv);

double String_Propagetor(double *z_p,double *z,int numcv,double **M,double *theta,double kappa,double gamma,double dt);

double String_Repara(double **z,double **z_p,int numreplica,int numcv);

double String_cTheta_FASYS(double *theta,double *crd);

double String_cM_FASYS(double *crd,double *mass,double **M, int numatom,int numcv);

#endif
