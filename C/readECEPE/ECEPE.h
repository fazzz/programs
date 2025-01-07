
#include <glib.h>

#define numatomtype 28
#define numatomtype2 50
#define numatomtype3 10
#define numatomtype4 4

struct ECEPE_dihe {
  double angle;
  int indexv1,indexv2;
  int ibnd1,ibnd2;

  int ifront;
  int ibchar1,ibchar2,ibchar3;

  int dpairs[4];

  double A;
  int NB,NS,IFTOR;
};

struct ECEPE_atom  {
  double refcoord[3];
  double charge;
  int nbtype;
  int kunit;
  int katom;
  int jatom;
  char name_atom[4];
  char name_res[3];

};

struct ECEPE_parms {
  int NUMATM;
  int NUMVAR;
  int NUMRES;
  int NUMINT;
  int NUMS;

  struct ECEPE_dihe *dihed;
  struct ECEPE_atom *atom;

};

struct pnb {
  double *Acff;
  double *Bcff;
};

struct ECEPE_pote {
  double p_t;
  double p_tors;
  double p_es;
  double p_nb;

};

struct ECEPE_force {
  double f_t[3];
  double n_tors;
  double f_es[3];
  double f_nb[3];

};


void read_ECEPE_parm(char *preofilename, char *bd8filename, struct ECEPE_parms *ECEPE_p, struct pnb *nb_p);

void read_ECEPE_parm_wtransindex(char *preofilename, char *bd8filename, struct ECEPE_parms *ECEPE_p, struct pnb *nb_p);

double calc_ff_ECEPE(struct ECEPE_pote *p, struct ECEPE_force *f, struct ECEPE_parms ECEPE_p, struct pnb nb_p, int *pairs, int numint);

double calc_TORS_ECEPE(double *t,double *n,double *A,int *ns,int *nb,double *dih, int *iftors,int numdih );

double calc_TORS_ECEPE2(double *t,double *n,double *A,int *ns,int *nb, int *iftors,int numdih, double *co, struct ECEPE_parms p );

double calc_NB_ECEPE(double *p_nb, double *f_nb,double *p_es,double *f_es,double *co,int *pairs,int numint,double *charge, int *nbtype,struct pnb nb_p);

int set_pairs_ECEPE(int *pairs,int numint,int *nbtype);

void read_ECEPE_coo (FILE *file, double *co, double *dihed, int numatom);

int make_int_pair_list(int **bp,int *numb,int numatom,int numbond, int **pair1_5, int **pair1_4, int *num14,int *num1_5);

int make_bd_pair_list(struct ECEPE_parms pa,int *bp,int *bp_f , int *numb);

void make_dihed_pairs_list( struct ECEPE_parms p, int **bpl, int *numb );

void make_dihed_pairs_list_v2( struct ECEPE_parms p, int **bpl, int *numb );

void sea_pair_by_name(char *name, int nb, int *bp,int ap[2], int apflag,GHashTable *bplist,GHashTable *bplist_c2,GHashTable *bplist_c3,GHashTable *bplist_c4,GHashTable *bplist_c5,struct ECEPE_parms pa);

int order_bd_pair_list(int *bp, int **bp_f, int numatom, int *numb;);

int lookup_ex_list(int **pair_ex_1_5,int numi, int numj, int *numex15);

double calc_NB_ECEPE_for_db(double *p_nb, double *f_nb,double *p_es,double *f_es,double *co,int **pairs_1_5,int **pairs_1_4,int *numnb,int *num14,int numatom,double *charge, int *nbtype,struct pnb nb_p);

double calc_ff_ECEPE_for_db(struct ECEPE_pote *p, struct ECEPE_force *f, struct ECEPE_parms ECEPE_p, struct pnb nb_p, int **pairs_1_5,int **pairs_1_4, int *num1_5, int *num1_4);

void read_ECEPE_detail_coo (FILE *file, double *co, int numatom);

double calc_TORS_for_check(double *co, int atom1, int atom2, int atom3, int atom4);

int read_ECEPE_detail_coo_cyc (FILE *file, double *co, double *dihed,double *ene);

double calc_TORS_ECEPE2_for_check(double *t,double *n,double *A,int *ns,int *nb, int *iftors,int numdih, double *co, struct ECEPE_parms p ,FILE *file, double *delta_dihed);

double calc_ff_ECEPE_for_db_cyc(struct ECEPE_pote *p, struct ECEPE_force *f, struct ECEPE_parms ECEPE_p, struct pnb nb_p, int **pairs_1_5,int **pairs_1_4, int *num1_5, int *num1_4, FILE *file, double *delta_dihed);
 
double calc_TORS_for_get_sabun(int numdih, double *co, struct ECEPE_parms p ,double *dihed_inspidas, double *delta_dihed);

double  calc_TORS_ECEPE2_for_check_out(int numdih, double *co,struct ECEPE_parms p ,FILE *file, double *delta_dihed);

double chang_index_dihed_pairs(int numdih, struct ECEPE_parms p );
