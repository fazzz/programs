
#define MAXNUNITERATION 20

#define AA 0
#define CA 1
#define HV 2

#define MD 0
#define AMBER 1

void bf_trajectry_ncL(int numatom, int numstep, 
		      double *mass,double *coord_ref,
		      double *rmsd_trj,char *trjname,char *traj_bfname,
		      int MODE,int flago, int IOMODE);

void ave_coord_ncbL(double *coord_ref ,char *trjname,int numatom,int numstep,int interval,int MODE, int IOMODE);


