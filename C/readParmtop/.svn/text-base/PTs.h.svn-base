
#define MNATOM 2000
#define MAXNTYPES 100
#define MAXNRES 400
#define MAXNUMBND 2000
#define MAXNUMANG 2000
#define MAXNNPTRA /*1000*/2000
#define MAXNSPM /*1000*/2000
#define MAXNBPER /*1000*/2000
#define MAXNGPER /*1000*/2000
#define MAXLESNTYPE 1000
#define MAXNATYP 1000
#define MAXNBONH 1000
#define MAXNBONA 1000
#define MAXNTHETH 2000
#define MAXNTHETA 2000
#define MAXNPHIH /*4000*/ 6000
#define MAXNPHIA 4000
#define MAXNEXT /*7000*/ 160000
#define MAXNPHB 1000
#define MAXNDPER 1000


struct AmberParm{
  char ITITLE[4];

  int NATOM;
  int NTYPES;
  int NBONH;
  int MBONA;
  int NTHETH;
  int MTHETA;
  int NPHIH;
  int MPHIA;
  int NHPARM;
  int NPARM;
  int NNB;
  int NEXT;
  int NRES;
  int NBONA;
  int NTHETA;
  int NPHIA;
  int NUMBND;
  int NUMANG;
  int NPTRA;
  int NATYP;
  int NPHB;
  int IFPERT;
  int NBPER;
  int NGPER;
  int NDPER;
  int MBPER;
  int MGPER;
  int MDPER;
  int IFBOX;
  int NMXPS;
  int IFCAP; 
  int NEXTRA;
  
  char *IGRAPH[4];
  double *CHRG;
  double *AMASS;
  int *IAC;
  int *NUMEX;
  int *ICO;
  char *LABERES[4];
  int *IPRES;
  double *RK;
  double *REQ;
  double *TK;
  double *TEQ;
  double *PK;
  double *PN;
  double *PHASE;
  double *SOLTY;
  double *CN1;
  double *CN2;
  int *BH[3];
  int *BA[3];
  int *TH[4];
  int *TA[4];
  int *PH[5];
  int *PA[5];
  int *NATEX;
  double *ASOL;
  double *BSOL;
  double *HBCUT;
  char *ISYMBL[4];
  char *ITREE[4];
  int *JOIN;
  int *IROTAT;
  int IPTRES;
  int NSPM;
  int NSPSOL;
  int *NSP;
  double BETA;
  double BOX[3];
  int NATCAP;
  double CUTCAP;
  double XCAP;
  double YCAP;
  double ZCAP;
  int BPER[MAXNBPER][2];
  int ICBPER[MAXNBPER*2];
  int TPER[MAXNGPER][3];
  int ICTPER[MAXNBPER*2];
  int PPER[MAXNDPER][4];
  int ICPPER[MAXNDPER*2];
  char IGRPER[MNATOM][4];
  char ISMPER[MNATOM][4];
  double ALMPER[MNATOM];
  double IAPER[MNATOM];
  int IACPER[MNATOM];
  double CGPER[MNATOM];
  double ATPOL[MNATOM];
  double ATPOL1[MNATOM];
  int NLES_NTYP;
  int LES_TYPE[MNATOM];
  double LES_FAC[MAXLESNTYPE*MAXLESNTYPE];
  double LES_CNUM[MNATOM];
  double LES_ID[MNATOM];


};

struct AmberParm AP;

int readParmtop(FILE *parmfile/*, struct AmberParm AP*/);
int readdihedpairs(int **atomdihedpairs, int *num);
int joinatomtores(int numatom, char LABERES[4]);
int writeParmtop(FILE *parmfile);
int sdvdWradii(double s);
int sd(double s);
int sd_woba(double s);
