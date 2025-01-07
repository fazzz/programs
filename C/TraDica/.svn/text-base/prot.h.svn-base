/////////////////////////////////////////
//                                     //
//          /MD_NA/src/ABA.h           //
//                                     //
/////////////////////////////////////////


// 定数

#define MAXNAC /*20*//*300*/600 /* 剛体数の上限*/
#define MAXA /*1000*/2000 /* 原子数の上限*/
#define MAXDOF /*300*/600 /* 自由度の上限*/
#define MAXNBH 4 /* 分岐数の上限*/
#define MAXNDAC /*10*/20 /* 剛体中の原子数の上限*/
#define MAXPDB 10000 /* 出力 PDB ファイル数の上限*/

// 変数
int IndexOfABICycle[MAXDOF];

double delta_dihed_time[MAXDOF][MAXNBH];
double dihed[MAXDOF][MAXNBH];

double dd_thetat[MAXDOF][MAXNBH];

double delta_matrix[4][4];

double Matrix[100][100];
double Hatrix[100];
double Patrix[100];

struct{
	int origin_atom_a;
	int terminal_atom_a[MAXNBH];
	int terminal;
	int num_clust;
	int num_atom_clust;
	int num_branch;
  int join;/*0410*/

	int nNumClutOfParent;
	int nNumClutOfChild[MAXNBH];

	int origin_xoord_a;
	int num_xoord_a;

	double xoord_clust[MAXNBH][MAXA][3];
	double mass_clust[MAXNAC];
	double sum_mass;
	double Inertia_clust[3][3];
	double InertiaMatrix[6][6];
	double Coriolis_acc_Mat[3][3];
	double momentum_clust;
	double Inertia_clust_total;
/////////////////////////////////////////////////////////////////////////////////////////
	double PsedoInertia[4][4];
	double PsedoTransMatrix[4][4];
	double qCOM[3];
/////////////////////////////////////////////////////////////////////////////////////////
	double dihedang[MAXNBH];
	double ddihedang[MAXNBH];
	double dddihedang[MAXNBH];
//	double old_dddihedang;
	double correct_dihedang[6];
	double predict_dihedang[6];
	double now_deltadihedang[MAXNBH];
	double Hing[MAXNBH][6];
	double TransMatrix[MAXNBH][6][6];
	double trans_A_to_CN[MAXNBH][3][3];

} clust[MAXDOF];





// 座標変換の関数
void trans_A_to_CN(int nNumClut);
void trans_CN_to_A(int nNumClut, int nNumClutOrigBranch);

int setJoin(int nNumClut, int joinflag);
