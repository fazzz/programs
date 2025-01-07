//////////////////////////////////////////
//                                      //
//          /MD_NA/src/force.h/         //
//                                      //
//////////////////////////////////////////

// íËêî
#define MAXDA 2000
#define MAXDIHED /*7000 (080510) */ 10000

int restflag;

// ïœêî
int kk;

//int atom_dihed_pair[MAXDIHED][6];
int *atom_dihed_pair/*[MAXDIHED][6]*/;

double V_dihed[MAXDIHED];
int n_dihed[MAXDIHED];
double theta_dihed[MAXDIHED];

// êßå¿óÕÇÃçÄ
double V_rest[MAXDIHED];
double theta_ref[MAXDIHED];

int dihed_rest[MAXDIHED][5];

int num_a_prot;
int NUM_A_PROT;

double q[MAXA][MAXA][3];
double len_q[MAXA][MAXA];
double F_L_J[MAXA][3];
double q_ele[MAXA];

double elect_pote[MAXA][MAXA]/*[100][100]*/;

double eata_i_j;
double ganmma_i_j;

int gnumnb;
int *gindexnb;
int gnum14;
int *gindex14;

int PEPCAACCMODE;
double fact;

// ä÷êî
void calc_force(int pflag);
void calc_force2(int pflag,
		 double *ele, double *ALJ, double *BLJ,
		 double *p_e, double *p_1_4_e, 
		 double *p_LJ,double *p_1_4_LJ,
		 double *f_e, double *f_1_4_e, 
		 double *f_LJ,double *f_1_4_LJ,
		 int numnb, int *indexnb,
		 int num14, int *index14, double *eig,double *eig_14, double *eig_dihed);

double Calc_L_J_PotentialandForce(void);
double Calc_L_J_P(int nNumClut,
	              int nNumClut2,
	              int i_c,
	              int ii_c);
void Calc_L_J_F(int nNumClut,
	            int i_c,
	            int j_c);
void Calc_L_J_F_1_4(int nNumClut,
	                int i_c,
	                int j_c);

double Calc_ele_sta_PotentialandForce(void);
double Calc_stat_elect_pote(int nNumClut,
	                        int nNumClut2,
	                        int i_c,
	                        int ii_c);
void Calc_stat_elect_f(int nNumClut,
	                   int nNumClut2,
	                   int i_c,
	                   int ii_c);
void Calc_stat_elect_f_1_4(int nNumClut,
	                       int nNumClut2,
	                       int i_c,
	                       int ii_c);

void Calc_dihed_Potential_Force(double *eig_dihed);
void Calc_dihed_Potential_Force_for_db(void);
double pick_dihed_one_clust(int nNumAom1,
	                      int nNumAom2,
	                      int nNumAom3,
	                      int nNumAom4,
	                      int nNumDihed);

int which_calc_non_bonding_pote(int nNumClut,
	                            int i_c,
	                            int j_a);
int which_calc_1_4_non_bonding_pote(int nNumClut,
	                                int i_c,
	                                int j_a);
void Calc_vector_atoms();
double power(double num,
	         int n);


double Calc_L_J_PotentialandForce2(double *ele, double *ALJ, double *BLJ,
				   double *p_e, double *p_1_4_e, 
				   double *p_LJ,double *p_1_4_LJ,
				   double *f_e, double *f_1_4_e, 
				   double *f_LJ,double *f_1_4_LJ,
				   int numnb, int *indexnb,
				   int num14, int *index14, double *cord, double *eig, double *eig_14);

void Calc_restraintforce(void);
void set_non_bonding_parameters(double *ele, double *ALJ, double *BLJ);
void set_non_bonding_index(void);

int inpindex[100];
int inpnum;
