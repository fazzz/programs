#ifndef INCLUDE_ABAb
#define INCLUDE_ABAb

struct clustdatab{
  int origin_atom_a;
  int *terminal_atom_a; // 2012-01-27
  int terminal;
  int num_clust;
  int num_atom_clust;
  int num_branch;
  int join;

  int nNumClutOfParent;
  int *nNumClutOfChild/*[4]*/;
  
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

typedef struct clustdatab CLTb;

struct ABAbdata{
  double corABIb[6][6],preABIb[6][6],KG[6],D;

  double corBF[6],preBF[6],nyu,eata;
};

typedef struct ABAbdata ABIb;

#define TERM 0

#include "ABAb_mainpass.h"
#include "ABAb_backpass.h"
#include "ABAb_prepass.h"
#include "ABAb_pick_data.h"
#include "ABAb_set_imat.h"
#include "ABAb_set_lref.h"
#include "ABAb_set_tmat.h"
#include "ABAb_set_trans.h"
#include "ABAb_set_frc.h"
#include "ABAb_update.h"
#include "ABAb_integ.h"
#include "ABAb_calcattfrc.h"
#include "ABAb_Inverse.h"
#include "ABAb_Inverse_backpass.h"
#include "ABAb_Inverse_mainpass.h"
#include "ABAb_gtree.h"
#include "ABAb_Nose-Hoover.h"
#include "ABAb_set_rst.h"
#include "ABAb_Nose-Hoover_new.h"
#include "ABAb_Nose-Hoover_new_mvV.h"
#include "ABAb_Nose-Hoover_chain.h"

#define NVT 1
#define NVE 0

#define k_B_kcm 1.98723e-3

double solverABAb(double *qacc,double *qvel,CLTb *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double q_NVT,double qvel_NVT,double *qacc_NVT,double s_NVT,double KE,double KEobj,int MODE);

double solverABAb_TermOn(double *qacc,double *qvel,CLTb *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double q_NVT,double qvel_NVT,double *qacc_NVT,double s_NVT,double *acc_Term,double *acc_Term2,double *vel_Term,double KE,double KEobj,int MODE);

double solverABAb_TermOn_NH(double *qacc,double *qvel,CLTb *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double s,double s_vel,double *s_acc,double tau2,double *acc_Term,double *acc_Term2,double *vel_Term,double Temp,double TempB);

double solverABAb_NH(double *qacc,double *qvel,CLTb *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double s,double s_vel,double *s_acc,double tau2,double Temp,double TempB);

double solverABAb_NH_new(double *qacc,double *qvel,CLTb *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double gzi,double *gzi_vel,double s,double *s_vel,double tau2,double Temp,double TempB);

double solverABAb_TermOn_NH_new(double *qacc,double *qvel,CLTb *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double gzi,double *gzi_vel,double s,double *s_vel,double tau2,double *acc_Term,double *acc_Term2,double *vel_Term,double Temp,double TempB);

double solverABAb_NH_new_mvV(double *qacc,double *qvel,CLTb *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double gzi);

double solverABAb_TermOn_NH_new_mvV(double *qacc,double *qvel,CLTb *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *q,double gzi,double *acc_Term,double *acc_Term2,double *vel_Term);

double solverABAb_NH_chain(double *qacc,double *qvel,CLTb *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *zeta_vel,double *zeta_acc,int M,int N,double KBT,double *Q_NH,double tau2,double Temp,double TempB);

double solverABAb_TermOn_NH_chain(double *qacc,double *qvel,CLTb *clt,double *Q,double *frc,double *crd,int numclt,int numatom,double *acc_Term,double *acc_Term2,double *vel_Term,double *zeta_vel,double *zeta_acc,int M,int N,double KBT,double *Q_NH,double tau2,double Temp,double TempB);

void ABAb_out_formated(FILE *outputfile,double pot,double KE,int i,double dt);

#endif

