
#define MAXNUNITERATION 20

#define AA 0
#define CA 1
#define HV 2

#define MD 0
#define AMBER 1


double bestfit(double *coord_ref,double *coord_tag,double *velo_tag,double *mass,int numatom,int flag);

void bf_trajectry(int numatom, int numstep, double *mass,double *coord_ref,double *rmsd_trj,char *trjname,char *velname, char *traj_bfname,char *velo_bfname,int flagv,int flaga,int flago, int flagc);

void ave_coord(double *coord_ref , FILE *inputfile, int numatom, int numstep, int flagc);

void bf_trajectry_nc(int numatom, int numstep, 
		     double *mass,double *coord_ref,
		     double *rmsd_trj,char *trjname,char *traj_bfname,
		     int MODE,int flago, int IOMODE);

void bf_trajectry_ncb(int numatom, int numstep, int interval,
		      double *mass,double *coord_ref,
		      double *rmsd_trj,char *trjname,char *traj_bfname,
		      int MODE,int flago, int IOMODE);

void ave_coord_nc(double *coord_ref ,char *trjname,int numatom,int numstep,int MODE, int IOMODE);

void ave_coord_ncb(double *coord_ref ,char *trjname,int numatom,int numstep,int interval,int MODE, int IOMODE);

double bestfit_nc(double *coord_ref,double *coord_tag,double *mass,int numatom);
