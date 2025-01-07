//////////////////////////////////////////
//                                      //
//          /MD_NA/src/MD.h/            //
//                                      //
//////////////////////////////////////////

// 定数
#define LIMITSTEPS 1000000

#define NVE 0

#define MAXDOF2 100

#define NVT 1

#define ON 1

#define OFF 0

#define LD 1

#define BD 2

#define MD 0

#define COO 1

#define CLU 2

#define TOP 3

#define SEQ 4

#define M_NUM 4

// 変数
double Cut_Off;
double deltat;
int TIME_LIMIT;
int NUM_IETRATION;

int MODE;
int MODEV;
int GearMODE;
int TermMoveMode;
int TermMoveMode2;
int MASSIVEOUT;

int restflag;
int dhstopflag;
int nbstopflag;
int veloutflag;

double Reatat;

int MAP;

int out_put_steps;
int out_put_steps_thomo;
int nNumStep;

char *InpfilCOORD;
char *InpfilCLUST;
char *InpfilSEQ;
char *InpfilTOP;
char *InpfilMDI;
char *InpfilCOORDREF;
char *InpfilDVELO;
char *InpfilEig;
char *OutfilTHMO;
char *OutfilCOORD;
char *OutfilVELO;
char *OutfilRESTAC;
char *OutfilRESTAV;
char *OutfilRESTAMFORM;
 
double old_deltadihedang[MAXDOF2];
double old_old_deltadihedang[MAXDOF2];

double acc_s_NVT;
double delta_s_NVT;

double GearsConstant[6];
double GearsConstant2[6];
double GearsConstant3[5];
double Telar_Matrix[6][6];
double Telar_Matrix2[5][5];



//double convertunitsys=1.660539e-27*1.0e-20/1.0e-24*2.3889e-4*6.022142e23;

// 関数
// velocity-verlet 法に関する関数
int velocity_verlet_step1(int n_step,int pflag,
			  double *ele, double *ALJ, double *BLJ,
			  double *p_e, double *p_1_4_e, 
			  double *p_LJ,double *p_1_4_LJ,
			  double *f_e, double *f_1_4_e, 
			  double *f_LJ,double *f_1_4_LJ,
			  int numnb, int *indexnb,
			  int num14, int *index14);

void velocity_verlet_step2(int pflag,			 
			   double *ele, double *ALJ, double *BLJ,
			   double *p_e, double *p_1_4_e, 
			   double *p_LJ,double *p_1_4_LJ,
			   double *f_e, double *f_1_4_e, 
			   double *f_LJ,double *f_1_4_LJ,
			   int numnb, int *indexnb,
			   int num14, int *index14,
			   double vel_Term[3], double *eig,double *eig_14, double *eig_dihed);

int verlet(int nNumStep);
// 2 面角加速度に関する関数
void calc_dd_theta_cycle(int pflag,int cflag,
			 double *ele, double *ALJ, double *BLJ,
			 double *p_e, double *p_1_4_e, 
			 double *p_LJ,double *p_1_4_LJ,
			 double *f_e, double *f_1_4_e, 
			 double *f_LJ,double *f_1_4_LJ,
			 int numnb, int *indexnb,
			 int num14, int *index14, double *eig,double *eig_14, double *eig_dihed);
void calc_sp_acc_cycle(int nNumClut);
void sub_calc_sp_acc_cycle(int nNumClut, int nNumCltminusone);
// 温度一定に関する関数
double velocity_scaling(double Energy_kinetic_o);
// マップ作成に関する関数
void CreatPsivsPhiMap(FILE *outputMap);


void getoption(int num, char *option);

void Gear(int nNumClut);
void Gear2(int nNumClut);

// NVT ensemble
int verlet_of_vartial_variable(int nNumStep);
// virtial variable
double xi_NVT;
// virtual variable parameter
double tau_NVT;
double q_NVT;
double s_NVT;
double dot_s_NVT;
//double dot_dot_s_NVT;
//double dot_dot_s_NVT_old;
double s_NVT_correct[6];
double s_NVT_predict[6];
double xi_NVT_correct[5];
double xi_NVT_predict[5];
double xi_dummy;
