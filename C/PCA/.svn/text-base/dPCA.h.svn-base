
#define MAXNUMDIHED 100

double *cov;
double *dihed_traj;
double *dihed_ave;
double *dihed_traj_trans;


void pick_dihed(/*double *dihed_traj,*/double *traj,int numdihed,int numatom, int atom_dihed_pair[MAXNUMDIHED][4], int time);

void calc_dPCA(/*double *dihed_traj,*/int time, int numdihed,
	       double covR[MAXNUMDIHED*MAXNUMDIHED], double w[MAXNUMDIHED]);

void calc_dcov(/*double *dihed_traj, */int time, int numdihed);

void proj_dPCA(/*double *dihed_traj,*/double *dihed_ave,int time, int numdihed,
	       double covR[MAXNUMDIHED*MAXNUMDIHED], double w[MAXNUMDIHED]);

