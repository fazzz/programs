
int get_path_ftrj(char *trjname,
		  char *pathnamebase,
		  double *crdref1, 
		  double *crdref2,
		  double criteria, 
		  double *path, 
		  int *numpoints,
		  int MODE, 
		  int interval);

int get_path_ftrj_bydha(double *trjname,
			char *pathnamebase,
			double *crdref1, 
			double *crdref2,
			double criteria,
			double *dha_spe,
			int *atomp,
			int numdha,
			double *path, 
			int *numpoints,
			int MODE, 
			int interval);

int execute_ini_path_wSBAAMC(double beta,double delta,
			     int ns,int interval,
			     char *ref1,char *ref2,char *crd,char *top,
			     char *dirbase,char *dirinipath,
			     char *inipathnc,
			     int numnode, char *name);


int anl_ini_path_SBAA(char *dirbase,char *dirinipath,char *inipathnc,
		      char *dirinipathanal, 
		      char *ene, 
		      char *dihed ,
		      char *rmsdf1 , char *rmsdf2, 
		      char *ref1, char *ref2,
		      int numnode, char *name);

int get_PATH_fSBAA(char *dirinipath,char *inipathnc,
		   double criteia,
		   char *cond, 
		   int numnode, char *name);

int modify_ini_path(char *dirinp,char *inipathbase,
		    int numpath,int numpoint,
		    char *dirout,char *modpsthbase,
		    int numnode, char *name);

int make_pdbfil_from_mpath(char *dirinp,char *modpsthbase,
			   int numpath, int numpoint,
			   char *dirout,char *pdbbase,
			   int numnode, char *name);

