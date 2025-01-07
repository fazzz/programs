#ifndef INCLUDE_ABA
#define INCLUDE_ABA

struct clustdata{
  int origin_atom_a;
  int terminal_atom_a[4];
  int terminal;
  int num_clust;
  int num_atom_clust;
  int num_branch;
  int join;

  int nNumClutOfParent;
  int nNumClutOfChild[4];
  
  double *xoord;
  double *mass;
  double sum_mass;
  double qCOM[3];
  double I[3][3];
  double IM[6][6];
  double trans_A_to_CN[3][3];

  double TM[6][6];

  double Spvel[6];
  double Spacc[6];
  double preSpacc[6];
  double Spfrc[6];
  double *frc_e;
  double *frc_LJ;
  double *frc_1_4_e;
  double *frc_1_4_LJ;
  double Coacc[6];
  double Cofrc[6];

  double F[6];
};

typedef struct clustdata CLT;

struct ABAdata{
  double corABI[6][6],preABI[6][6],KG[6],D;

  double corBF[6],preBF[6],nyu,eata;
};

typedef struct ABAdata ABI;

#define TERM 0

#include "ABA_mainpass.h"
#include "ABA_backpass.h"
#include "ABA_prepass.h"
#include "ABA_pick_data.h"
#include "ABA_set_imat.h"
#include "ABA_set_lref.h"
#include "ABA_set_tmat.h"
#include "ABA_set_trans.h"
#include "ABA_set_frc.h"
#include "ABA_update.h"
#include "ABA_integ.h"
#include "ABA_calcattfrc.h"
#include "ABA_Inverse.h"
#include "ABA_Inverse_backpass.h"
#include "ABA_Inverse_mainpass.h"
#include "ABA_gtree.h"
#include "ABA_Nose-Hoover.h"
#include "ABA_set_rst.h"
#include "ABA_Nose-Hoover_new.h"

#define NVT 1
#define NVE 0

#define k_B_kcm 1.98723e-3

double solverABA(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double q_NVT,double qvel_NVT,double *qacc_NVT,double s_NVT,double KE,double KEobj,int MODE);

double solverABA_TermOn(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double q_NVT,double qvel_NVT,double *qacc_NVT,double s_NVT,double *acc_Term,double *acc_Term2,double *vel_Term,double KE,double KEobj,int MODE);

double solverABA_TermOn_NH(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double s,double s_vel,double *s_acc,double tau2,double *acc_Term,double *acc_Term2,double *vel_Term,double Temp,double TempB);

double solverABA_NH(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double s,double s_vel,double *s_acc,double tau2,double Temp,double TempB);

double solverABA_NH_new(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double gzi,double *gzi_vel,double s,double *s_vel,double tau2,double Temp,double TempB);

double solverABA_TermOn_NH_new(double *qacc,double *qvel,CLT *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double gzi,double *gzi_vel,double s,double *s_vel,double tau2,double *acc_Term,double *acc_Term2,double *vel_Term,double Temp,double TempB);

void ABA_out_formated(FILE *outputfile,double pot,double KE,int i,double dt);

#endif
