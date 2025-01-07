
#define ON 1
#define OFF 0

int Create_inpcrd_from_trj(char *outputfilenamebase,FILE *inputfile,int numsnap, int numinpcrd, int numatom,int *num, int amberflag);

int Create_inpcrd_from_trj_Amber(char *outputfilenamebase,FILE *inputfile,int numsnap, int numinpcrd, int numatom,int *num, int amberflag);

int Create_inpcrd_from_trj_with_ene_check(char *outputfilenamebase, char *inputfilename,int numsnap, int numinpcrd, int numatom,int *num, int amberflag, double *ene, double checkv);
