/////////////////////////////////////////
//                                     //
//          /MD_NA/src/ABA.h           //
//                                     //
/////////////////////////////////////////


#define MAXNAC 100
#define MAXDOF 1000
#define MAXNBH 4 
#define MAXNDAC 20
#define MAXPDB 10000

int *IndexOfABICycle;

double **delta_dihed_time;
double **dihed;

double **dd_thetat;

double delta_matrix[4][4];

double Matrix[100][100];
double Hatrix[100];
double Patrix[100];

double delta_Term[6],vel_Term[6],acc_Term[6],acc_Term2[6],acc_Term3[6];
double vel_q_Term[4],Rot_Term[3][3];
double Euler[3],dotEuler[3];

// Articulated-Body Inertia
struct Articulated_Body_Inertia
{
	double PredictABI[6][6];
	double CorrectABI[6][6];
	double ABJointI;
	double Kalman_Gain[6];
	double Kalman_Gain_Transpose[6];
	int hingmat;
}*ABI;

// Bias Force
struct Bias_force
{
	double Predictzzz[6];
	double Correctzzz[6];
	double eata;
	double nyu;
	double hingforce;
}*zzz;

struct num_ATOM_belong_dihed
{
	int num_ATOM_N;
	int num_ATOM_N_1;
	double dihedang;
};

struct non_bonding_pairs
{
  int **not_interacting;
};

struct o_f_non_bonding_pairs
{
  int **o_f_not_interacting;
};

struct force_parm
{
  double *e_f;

  double *V_dihed;
  double *theta_dihed;
  int *n_dihed;
  int num_dihed_clust;
  int num_dihed_clust_1;
  int num_dihed_clust_2;
  int num_ATOM_N[10];
  int num_ATOM_N_1[10];

  struct num_ATOM_belong_dihed n_A_B_dihed[10];
};

// spatial force
struct spatial_force
{
  double N_clust[3];
  double f_clust[3];
};

// ���Τ���
struct force
{
  double bond;
  double angle;

  double f_dihed;
  double f_rest;
  double **f_elesta;
  double **f_L_J;
  double **f_1_4_elesta;
  double **f_1_4_L_J;

  double **f_atoms;
  double f_Brownian[6];

  struct spatial_force sp_f_clust[MAXNBH];
};

struct potential
{
  double bond;
  double angle;

  double p_dihedt;
  double p_elestat;
  double p_L_Jt;
  double p_1_4_elestat;
  double p_1_4_L_Jt;
  
  double p_restt;
  
  double p_total;

  double *p_elesta;
  double *p_L_J;
  double *p_1_4_elesta;
  double *p_1_4_L_J;
  double *p_at;
  double *p_dihedc;
  double *p_rest;

} potential_pro;

struct clustdata{
  int origin_atom_a;
  int terminal_atom_a[MAXNBH];
  int terminal;
  int num_clust;
  int num_atom_clust;
  int num_branch;
  int join;

  int nNumClutOfParent;
  int nNumClutOfChild[MAXNBH];
  
  int origin_xoord_a;
  int num_xoord_a;

  double **xoord_clust/*[MAXNBH][MAXA][3]*/;
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
  
  double sp_velo[6];
  double Coriolis_acc[6];
  double Coriolis_b[6];
  double sp_acc[6];
  double predict_alpha[6];
  double predict_velo[6];

  struct force_parm f_p_clust;
  struct force f_c;
  struct non_bonding_pairs pairs;
  struct o_f_non_bonding_pairs o_f_pairs;
} *clust;

void calc_ABA_cycle(int nNumClut);
void sub_calc_initial_ABI(int nNumClut);
void sub_calc_ABI_cycle(int nNumCluto,
			int nNumClutplusone[10],
                        int num_branch);
void sub_sub_calc_ABI_cycle(int nNumCluto, int nNumClutplusone);
int calcPredictABI(int nNumCluto);
void calcKalmanGain(int nNumCluto);
void CalcD(int nNumCluto);
void sub_CalcCorrectABI(int nNumCluto, int nNumClutplusone);
void CalcCorrectABI(int nNumCluto, 
		    int nNumClutplusone[10], 
		    int num_branch);

void calc_Bias_force_cycle(int nNumClut);
void sub_calc_Bias_force_initial(int nNumClut);
void sub_calc_Bias_force_cycle(int nNumCluto, 
			       int nNumClutplusone[10],
			       int num_branch);
void sub_sub_calc_Bias_force_cycle(int nNumCluto, int nNumClutplusone);
void calcPredictzzz(int nNumCluto);
void calceata(int nNumCluto);
void calcnyu(int nNumCluto);
void CalcCorrectzzz(int nNumCluto, 
	                int nNumClutplusone[10], 
	                int num_branch);

void pick_dihed(int nNumClut);

void set_atom_velo(void);
void sub_set_atom_velo(int nNumClut, int nNumBod);
void set_coriolis_acc(int nNumClut);
void sub_set_coriolis_acc(int k_o, int kk, int nb);
void set_sp_velo(int nNumClut, int nNumClutOrigBranch);
void sub_set_sp_velo(int nNumClut, int nNumClutminusone);
void set_coriolis_force(int nNumClut);
void sub_set_coriolis_force(int nNumClut, int nNumClutminousone);
void set_predict_acc(int nNumClt, int nNumClutOrigBranch);
void sub_set_predict_acc(int nNumClt, int nNumCltminusone);

void trans_A_to_CN(int nNumClut);
void sub_trans_A_to_CN(int nNumCltTar, int nNumCltCoo,
                       int nNumBod, int nNumAtom);
void sub_trans_A_to_CN_Initial(void);
void trans_CN_to_A(int nNumClut, int nNumClutOrigBranch, double q_Term[4],double dq_Term[6]);
void sub_trans_CN_to_A(int nNumCltTar, int nNumAtom, int nNumAtomLoca);

void output();
void output_file_coo_pro(FILE *output_c);
void output_file_vel_pro(FILE *output_v);
void output_file_dihed_vel_pro(FILE *output_v);
void output_file_initial_pdb(void);
void output_file_pdb(FILE *output_c, int nNumStep);
void output_restat(FILE *output_rest_1, 
                   FILE *output_rest_2, int nNumStep);

void Initialize_Variables(void);

void pick_dihed2(int nNumClut);



void calc_Matrix_cycle(void);
double sub_calc_Matrix_cycle(int nNumClut1, int nNUmClut2);
void sub_sub_calc_Matrix_cycle(int nNumClutI,
			       int nNumClutJ,
			       int nNumClutK);
void calc_dot_Pseduo_TransMatrix(int nNumClutI,
				 int nNumClutK,
				 double mat[4][4]);
void calc_Hatrix_cycle(void);
double sub_calc_Hatrix_cycle(int nNumClutI, int nNumClutJ, int nNumClutM);
double sub_sub_calc_Hatrix_cycle(int nNumClutI,
				 int nNumClutJ,
				 int nNumClutM,
				 int nNumClutK);

void calc_dot_dot_Pseduo_TransMatrix(int nNumClutI,
				     int nNumClutJ,
				     int nNumClutK,
				     double mat[4][4]);

void Jordan(int numOfRow, double Mat[100][101], double vect[100]);

double *old_dddihedang;

int setJoin(int nNumClut, int joinflag);
