//////////////////////////////////////////
//                                      //
//       /MD_NA/src/physics.h/          //
//                                      //
//////////////////////////////////////////

/* physical_constants*/
#define k_B_J 1.38065e-23
#define AB_num 6.022142e23
#define eata_o 8.8541878e-12
#define h_plnk 6.626069e-34
#define visco_wat 8.94e-4
#define kcaltoJ /*6.95110e-21*/2.3889e-4
#define Jtokcal 1.43862e21
#define mutokg 1.660539e-27
#define stofs 1.00e-15
#define stons 1.00e-12
#define stops 1.00e-9
#define mtoA 1.00e-10

#define k_B_kcm 1.98723e-3 //1.3807e-23/(4.184*1000.0)*6.022e23


/* thermodynamical properties*/
double T_Kelvin;
double Energy_Internal;
double Energy_kinetic;
double Energy_kinetic_o2;
double Energy_potential;
double V_system;

double Energy_kinetic8;

double Entorpy_system;

double Energy_kinetic_stac;
double Energy_potential_stac;
double Energy_Internal_stac;

// åªç›ÇÃâ∑ìx
double T_Kelvin_Now;
double T_Kelvin_Now2;
double T_Kelvin_Initial;

int DOFOFPROT;

/* transformation constants*/

/* functions*/
void set_initial_velo(void);
void set_initial_velode(void); // 0911
double analthermo_dyn_properties(FILE *output, double vel_Term[3]);
void CalcTotal_Energy(double vel_Term[3]);
void CalcT(double vel_Term[3]);
void calc_velo2(void);
