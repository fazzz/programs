//////////////////////////////////////////
//                                      //
//          /MD_NA/src/gener.h          //
//                                      //
//////////////////////////////////////////

// 定数
#define MAXRES 400 /* 残基数の上限*/
#define MAXPEPTIDE /*100 (080510)*/1000 /* ペプチド結合の上限*/
#define MAXA  /*2000 (080510)*/3000 /* 原子数の上限*/
#define MAX_2 /*500*//*600 (080510)*/1000 /* 自由度の上限*/
#define MAXA_RES 40 /* 残基中の原子数の上限*/
#define MAXNAME 40
#define MAXCHI 20
#define MAXL 10
#define KATS 10
#define PI 3.141592653589793238
#define AMBERMODE 1

#define TERMINAL 0
#define NOTTERMINAL 1
#define BRANCH 2

/* 変*/

int n;
int *npre/*[MAXRES]*/;
double ident[6][6];

int nNumAtomPeptide_N[MAXPEPTIDE];
int nNumAtomPeptide_C[MAXPEPTIDE];

struct L_J_parmat{
  double *A/*[10000]*/;
  double *B/*[10000]*/;

  int *atomtype/*[1000]*/;
  int **atomtypeIndex/*[1000][1000]*/;
  //  double A[10000];
  //  double B[10000];

  //  int atomtype[1000];
  //  int atomtypeIndex[1000][1000];


};

// タンパク質のデータ
struct protein{
	char *name_prot;
	int num_atom;
	int nNumPeptide;
  int *nNumClutPeptide/*[MAXPEPTIDE]*/;
	int inumrs;
	int num_clust;
	int DOF;
	int nNumDihedALL;
	int nNumDihedType;

        int nNumDihed_rest;

  double **coord/*[MAXA][3]*/;
	/*************************************/
        /* double old_coord[MAXA][3];	     */
        /*************************************/
  double **velo/*[MAXA][3]*/;
	/************************************/
        /* double old_velo[MAXA][3];	    */
	/* double acc[MAXA][3];		    */
        /************************************/
  double qCOM[3];
	double qCOM_Old[3];
	double veloCOM[3];
	double sumMass;
  int *name_atom/*[MAXA]*/;
	struct L_J_parmat L_J_parm;
}prot;

// 残基のデータ
struct res{
	int num_atom;
	int res_ID;
}residue[10];

double mq[3];

// 関数
// データの取得を行う関数
int pick_data(int pflag);

// 種種の設定を行う関数
void set(void);
//void set_Inertia_clust(int nNumClut);
void set_Inertia_clust(int nNumClut,int  num_atom_clust_total);
double Sum_Mass(int nNumClut);
void InertiaMatrix(int nNumClut);

void set_trans_Matrix(int nNumClt, int nNumClutOrigBranch);
void sub_set_trans_Matrix(int nNumCltplusone, int nNumClt);
void set_pseduo_trans_matrix(int nNumClut);
