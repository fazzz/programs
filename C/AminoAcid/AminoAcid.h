
#ifndef INCLUDE_AA
#define INCLUDE_AA

#define NHASH 43
#define MAXNRESIDUE 100
#define NRESIDUE 100     

struct AminAcidData{
  int numatom;
  char **atomname;

  int numdihed;  
  int numkai;
  char ***atomnamepair_kai;

  int numatomHead;
  int numatomTail;
};

typedef struct AminAcidData AAdata ;

struct AminAcidData dataofAAontable[25];

int NumResInStdResData;


struct AminAcidDatatable {
  char nameofthisresidue[20];
  
  AAdata AAdataontable;
  struct AminAcidDatatable *next;
};

typedef struct AminAcidDatatable AADatatable;

AAdata *table[NHASH];

int readAADataBase(char *AAdatafilename);
AADatatable *LURAAData(char name[10],int create,AAdata dataofthisAminoAcid);
unsigned int hash(char * str);

#endif
