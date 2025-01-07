

#ifndef INCLUDE_REMD_TAM2_2
#define INCLUDE_REMD_TAM2_2

#define AA1INPF 0
#define CG1INPF 1
#define CG2INPF 2
#define AA1KZ 3
#define CG1KZ 4
#define CG2KZ 5

double **MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_Amber_2PROTEINS2008_2013_08_31(int myrank,int num_procs,
									     int tag, MPI_Status* status,
									     int numRE, int numEX, double *KZAA, double **KZCG,
									     struct AADataforREMD_Amber AAData,
									     struct CGDataforREMD_PROTEINS2008 CGData[2],
									     struct TACCMDataforREMD_A_2P2008 ZData,
									     struct AACGCommonDataforREMD_A_P2008 CData,
									     double T0AA,double T0CG[2], double T0Z, 
									     int numstep, int interval, 
									     double dt,double dt2,
									     double wdt2[3],double wdt4[3], int nc,
									     double UNITT, double k_B, double tau, double pi,
									     FILE* logfile, 
									     double fact_b,double fact_a,double fact_t,double fact_NC,double fact_NNC);


#endif
