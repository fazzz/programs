
#define MAXNATOM 1000
#define MAXNTYPES 100
#define MAXNRES 100
#define MAXNUMBND 1000
#define MAXNUMANG 1000
#define MAXNNPTRA 1000
#define MAXNSPM 100
#define MAXNBPER 100
#define MAXNGPER 100
#define MAXLESNTYPE 100
#define MAXNATYP 100
#define MAXNBONH 100
#define MAXNBONA 100
#define MAXNTHETH 100
#define MAXNTHETA 100
#define MAXNPHIH 100
#define MAXNPHIA 100
#define MAXNEXT 100
#define MAXNPHB 100
#define MAXNDPER 100


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
  
  char IGRAPH[MAXNATOM][4];

  double CHRG[MAXNATOM];

  double AMASS[MAXNATOM];

  int IAC[MAXNATOM];

  int NUMEX[MAXNATOM];

  int ICO[MAXNTYPES*MAXNTYPES];

  char LABERES[MAXNRES][4];

  int IPRES[MAXNRES];

  double RK[MAXNUMBND];

  double REQ[MAXNUMBND];

  double TK[MAXNUMANG];

  double TEQ[MAXNUMANG];

  double PK[MAXNNPTRA];

  double PN[MAXNNPTRA];

  double PHASE[MAXNNPTRA];

  double SOLTY[MAXNATYP];

  double CN1[MAXNTYPES*(MAXNTYPES+1)/2];

  double CN2[MAXNTYPES*(MAXNTYPES+1)/2];

  int BH[MAXNBONH][3];

  int BA[MAXNBONA][3];

  int TH[MAXNTHETH][4];

  int TA[MAXNTHETA][4];

  int PH[MAXNPHIH][4];

  int PA[MAXNPHIA][4];

  int NATEX[MAXNEXT];

  double ASOL[MAXNPHB];

  double BSOL[MAXNPHB];

  double HBCUT[MAXNPHB];

  char ISYMBL[MAXNATOM][4];

  char ITREE[MAXNATOM][4];

  int JOIN[MAXNATOM];

  int IROTAT[MAXNATOM];

  int IPTRES;

  int NSPM;

  int NSPSOL;

  int NSP[MAXNSPM];

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

  char IGRPER[MAXNATOM][4];

  char ISMPER[MAXNATOM][4];

  double ALMPER[MAXNATOM];

  double IAPER[MAXNATOM];

  int IACPER[MAXNATOM];

  double CGPER[MAXNATOM];

  double ATPOL[MAXNATOM];

  double ATPOL1[MAXNATOM];

  int NLES_NTYP;

  int LES_TYPE[MAXNATOM];

  double LES_FAC[MAXLESNTYPE*MAXLESNTYPE];

  double LES_CNUM[MAXNATOM];

  double LES_ID[MAXNATOM];


};

struct AmberParm AP;

int readParmtop(FILE *parmfile/*, struct AmberParm AP*/);

