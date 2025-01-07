
#define ON 1
#define OFF 0

#define AA 0
#define CA 1
#define HV 2

typedef struct PDBatom PDBA;

struct PDBatom {
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

typedef struct PDBform PDBF;

struct PDBform {
  int numatom;
  PDBA *PDBa;
};

int writePDB(FILE *pdbfile,PDBF PDB);
int writPDB_wopt(FILE *pdbfile,PDBF PDB, int MODE);
int readPDB(FILE *pdbfile,PDBF PDB,int numatom);
int copyPDBform(PDBF PDB1,PDBF PDB2);
int readPDBatomnum(FILE *pdbfile,int numatom);
int readPDBdatafromParmtop(PDBF PDB);
