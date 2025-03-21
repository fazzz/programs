#ifndef INCLUDE_PTL
#define INCLUDE_PTL

#define MNATOM /*7000*/9000
#define MAXNTYPES 100
#define MAXNRES /*400*/600
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


struct AmberParmL{
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
  
  char IGRAPH[MNATOM][4];

  double *CHRG/*[MNATOM]*/;
  double *AMASS/*[MNATOM]*/;
  int *IAC/*[MNATOM]*/;
  int *NUMEX/*[MNATOM]*/;
  int *ICO/*[MAXNTYPES*MAXNTYPES]*/;
  char LABERES[MAXNRES][4];
  int *IPRES/*[MAXNRES]*/;
  double *RK/*[MAXNUMBND]*/;
  double *REQ/*[MAXNUMBND]*/;
  double *TK/*[MAXNUMANG]*/;
  double *TEQ/*[MAXNUMANG]*/;
  double *PK/*[MAXNNPTRA]*/;
  double *PN/*[MAXNNPTRA]*/;
  double *PHASE/*[MAXNNPTRA]*/;
  double *SOLTY/*[MAXNATYP]*/;
  double *CN1/*[MAXNTYPES*(MAXNTYPES+1)/2]*/;
  double *CN2/*[MAXNTYPES*(MAXNTYPES+1)/2]*/;
  int **BH/*[MAXNBONH][3]*/;
  int **BA/*[MAXNBONA][3]*/;
  int **TH/*[MAXNTHETH][4]*/;
  int **TA/*[MAXNTHETA][4]*/;
  int **PH/*[MAXNPHIH][5]*/;
  int **PA/*[MAXNPHIA][5]*/;
  int *NATEX/*[MAXNEXT]*/;
  double *ASOL/*[MAXNPHB]*/;
  double *BSOL/*[MAXNPHB]*/;
  double *HBCUT/*[MAXNPHB]*/;
  char ISYMBL[MNATOM][4];
  char ITREE[MNATOM][4];
  int *JOIN/*[MNATOM]*/;
  int *IROTAT/*[MNATOM]*/;
  int IPTRES;
  int NSPM;
  int NSPSOL;
  int *NSP/*[MAXNSPM]*/;
  double BETA;
  double BOX[3];
  int NATCAP;
  double CUTCAP;
  double XCAP;
  double YCAP;
  double ZCAP;
  int **BPER/*[MAXNBPER][2]*/;
  int *ICBPER/*[MAXNBPER*2]*/;
  int **TPER/*[MAXNGPER][3]*/;
  int *ICTPER/*[MAXNBPER*2]*/;
  int **PPER/*[MAXNDPER][4]*/;
  int *ICPPER/*[MAXNDPER*2]*/;
  char **IGRPER/*[MAXNATOM][4]*/;
  char **ISMPER/*[MNATOM][4]*/;
  double *ALMPER/*[MNATOM]*/;
  double *IAPER/*[MNATOM]*/;
  int *IACPER/*[MNATOM]*/;
  double *CGPER/*[MNATOM]*/;
  double *ATPOL/*[MNATOM]*/;
  double *ATPOL1/*[MNATOM]*/;
  int NLES_NTYP;
  int *LES_TYPE/*[MNATOM]*/;
  double *LES_FAC/*[MAXLESNTYPE*MAXLESNTYPE]*/;
  double *LES_CNUM/*[MNATOM]*/;
  double *LES_ID/*[MNATOM]*/;

};

struct AmberParmL AP;

int readParmtopL(FILE *parmfile);
int readdihedpairsL(int **atomdihedpairs, int *num);
int PTL_joinatomtores(int numatom, char LABERES[4]);

int PTL_res_ca(int numres);

//int PTL_iniatomnumofres(int numres);
int PTL_resnum(int numatom);

int PTL_resnum2(int numatom);

int PTL_canum_fromresnum(int numres);

int PTL_which_include(int numres,int *listres,int numlistres);

#endif
