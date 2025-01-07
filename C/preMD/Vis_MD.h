/////////////////////////////////////////
//                                     //
//         /Vis_MD/src/DataBase.h      //
//                                     //
/////////////////////////////////////////

#define MAXNATOM 4000/*9000*/
#define MAXNCLUST /*2000*/3000
#define MAXNDH 1000    
#define MAXNTD 400       
#define MAXNRESIDUE 1000 
#define MAXNBH 4        
#define NRESIDUE 100     
#define NHASH /*26*//*29*//*39*/43
#define MMAXNATOM  400
#define MAXNATOMINCLUT 100
#define MAXNATOMINP 4000
#define MAXNCLUSTINP 3000
#define MAXNDHINP /*1000*//*4000*/10000

#define PI 3.141592653589793238

char *InpfilCOORD;
char *InpfilCLUST;
char *InpfilSEQ;
char *InpfilTOP;

int ADDFLAG;

typedef struct atomtype atomtype;
struct atomtype {
  int atomnum;
  double gamma;
  double eata;
  double e_elect;
};

typedef struct atomdata atomdata;
struct atomdata {
  struct atomtype atomID;
  double coord[3];
  
  int bond_interacting[10];
  int num_bond_interacting;

  int not_interacting[50];
  int num_not_interacting;
  int o_f_interacting[50];
  int num_o_f_interacting;
};


typedef struct dihedang dihedang;
struct dihedang{
  int atomnum[6];
};

typedef struct dihedang_parm dihedang_parm;
struct dihedang_parm{
  double V_2;
  double theta;
  int num;
};

typedef struct dataofdiheds dataofdiheds;
struct dataofdiheds{
  int atomnum[4];
  
  struct dihedang_parm dihedangp;
};

struct dataofdiheds dataofdihed[20];

typedef struct clustdata clustdata;
struct clustdata{
  int clustID;
  int origin;
  int term[MAXNBH];

  int termflag;
  int nNumAtom;
  int nNumBranch;
  int nDegOfFoo;
  int nNumChildClust[MAXNBH];
  int nNumParentClust;

  int nNumDihed;
  int nNumDiehd_f[MAXNBH];
  int nNumDihed_b;
  int nNumAtom_f;
  int nNumAtom_b;
  
};

typedef struct residuedata residuedata;
struct residuedata{
  int nNumAtom;
  int nNumClut;
  int nNumDihed;
  
  struct atomdata atm[MAXNATOM];
  struct clustdata clt[MAXNCLUST];
  struct dihedang dihedang[MAXNDH];

  int HeadAtom;
  int TailAtom;
  int HeadClut;
  int TailClut;
};

struct residuedata dataofresidueontable[/*26*/40];

int NumResInStdResData;

typedef struct elementoftable elementoftable;

struct elementoftable {
  char nameofthisresidue[20];
  
  residuedata dataofthisresdiue;
  elementoftable *next;
};

elementoftable *table[NHASH];

double xoord[MAXNATOM][3];

double MatTransAbsoToLocal[3][3];

struct protein{
  char *Sequence[MAXNRESIDUE];
  int numsequence[MAXNRESIDUE];
  char *NameProt;

  int totalatomtype;
  int nNumAtom;
  int nNumAtomType;
  int nNumClut;
  int nNumResidue;
  int nNumDihedAng;
  int nNumDihedType;
  double atmtypeIndex[100];
  double parmA[100];
  double parmB[100];

  struct atomdata atm[MAXNATOMINP];
  struct atomtype atmtype[MAXNATOMINP];
  struct clustdata clt[MAXNCLUSTINP];
  struct dihedang dihedang[MAXNDHINP];
  struct dihedang_parm dihedang_parm[MAXNDHINP];

}prot;

int PickData(void);

int CreateData(void);

void TransAbsoToLocal(int nNumAtom1, int nNumAtom2, int nSumNumAtom);
void TransLocalToAbso(int nNumAtom1, int nNumAtom2, int nSumNumAtom, 
                      double delta_dihed);

void CreateInPut(int TotalAtomType);
void CreateCoordFile(FILE *input);
void CreateSeqFile(FILE *input);
void CreateClustFile(FILE *input);
void CreateTopFile(FILE *input, int TotalAtomType);

elementoftable *LURData(char nameo[10],int create,residuedata dataofthisresdiue);

unsigned int hash(char * str);

int MatchResName(char *nameofthisresidue);

void *emalloc(size_t n);
FILE *efopen(char *filename,char *flag);

void usage (char *progname);
///////////////////////////////////////////////////////////////////////////////


