
#define ON 1
#define OFF 0

#define AA 0
#define CA 1
#define HV 2

typedef struct PDBLatom PDBLA;

struct PDBLatom {
  int HETEROflag;
  int serial;
  char name[4];
  char altLOC;
  char resname[3];
  char ChainID;
  int resSeq;
  char iCode;
  double coord[3];
  double occupancy;
  double tempfact;
  char element[2];
  char charge[2];
};

typedef struct PDBLform PDBLF;

struct PDBLform {
  int numatom;
  PDBLA *PDBLa;
};

int writePDBL(FILE *pdbfile,PDBLF PDBL);
int writPDBL_wopt(FILE *pdbfile,PDBLF PDBL, int MODE);
int writPDBL_wopt_series(FILE *pdbfile,PDBLF PDBL, int MODE);
int readPDBL(FILE *pdbfile,PDBLF PDBL,int numatom);
int copyPDBLform(PDBLF PDBL1,PDBLF PDBL2);
int readPDBLatomnum(FILE *pdbfile,int *numatom);
int readPDBLdatafromParmtop(PDBLF PDBL);
int addCAPL(PDBLF PDBL1,PDBLF PDBL2);
//int findrespoint(PDBF PDB);

