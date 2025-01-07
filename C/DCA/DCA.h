#ifndef INCLUDE_DCA
#define INCLUDE_DCA

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
  double I[3][3];
  double IM[6][6];
  double trans_A_to_CN[3][3];

  double TransMatrix[6][6];

  double Spvel[6];
  double Spfrc[6];
  double *frc;
  double *Coacc;
  double *Cofrc;
};

typedef struct clustdata CLT;

struct DCAdata{
  double *P1,*P2,*P12,*P21;
  double *f1,*f2;
  double *b1,*b2;
  double *gamma,*W,*beta,*V;
};

typedef struct DCAdata DCA;

typedef struct assembletreedata AST;

struct assembletreedata {
  int leafflag;
  int right;
  int left;
  int refright;
  int refleft;
  int num;
};

#include "DCA_mainpass.h"
#include "DCA_backpass.h"
//#include "DCA_prepass.h"
#include "DCA_pick_data.h"
#include "DCA_set_imat.h"
#include "DCA_set_lref.h"
#include "DCA_set_tmat.h"
#include "DCAs_trans.h"
#include "DCAs_setfrc.h"

#endif
