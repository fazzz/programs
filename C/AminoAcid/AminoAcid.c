#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "EF.h"
#include "AminoAcid.h"

#define MULTIPLIER 31

int readAADataBase(char *AAdatafilename) {
  int i,j;
  int nNumRes;
  char	nameofthisresidue[MAXNRESIDUE][4];
  int	flag;

  char *line;
  size_t len=0;

  FILE *AAdatafile;

  AAdatafile=efopen(AAdatafilename,"r");
  
  for (nNumRes=0 ; nNumRes < NRESIDUE ; ++nNumRes){
    fscanf(AAdatafile, "%s", &(nameofthisresidue[nNumRes]));
    nameofthisresidue[nNumRes][3] = '\0';

    fscanf(AAdatafile, "%d", &(dataofAAontable[nNumRes].numatom));
    dataofAAontable[nNumRes].atomname=(char **)gcemalloc(sizeof(char *)*dataofAAontable[nNumRes].numatom);
    for (i=0;i<dataofAAontable[nNumRes].numatom;++i) {
      dataofAAontable[nNumRes].atomname[i]=(char *)gcemalloc(sizeof(char)*4);
    }

    for (i=0;i<dataofAAontable[nNumRes].numatom;++i) {
      fscanf(AAdatafile, "%s", dataofAAontable[nNumRes].atomname[i]);
    }

    fscanf(AAdatafile, "%d", &(dataofAAontable[nNumRes].numdihed));
    fscanf(AAdatafile, "%d", &(dataofAAontable[nNumRes].numkai));
    dataofAAontable[nNumRes].atomnamepair_kai=(char ***)gcemalloc(sizeof(char **)*dataofAAontable[nNumRes].numatom);
    for (i=0;i<dataofAAontable[nNumRes].numatom;++i) {
      dataofAAontable[nNumRes].atomnamepair_kai[i]=(char **)gcemalloc(sizeof(char *)*5);
      for (j=0;j<5;++j)
	dataofAAontable[nNumRes].atomnamepair_kai[i][j]=(char *)gcemalloc(sizeof(char )*4);
    }

    for (i=0;i<dataofAAontable[nNumRes].numkai;++i) {
      for (j=0;j<5;++j) {
	fscanf(AAdatafile, "%s", dataofAAontable[nNumRes].atomnamepair_kai[i][j]);
      }
    }

    //    fscanf(AAdatafile, "%d", &(dataofAAontable[nNumRes].numatomHead));
    //    fscanf(AAdatafile, "%d", &(dataofAAontable[nNumRes].numatomTail));

    LURAAData(nameofthisresidue[nNumRes],0,dataofAAontable[nNumRes]);
    getline(&line,&len,AAdatafile);
    getline(&line,&len,AAdatafile);

    if (++NumResInStdResData > NRESIDUE) {
      printf("too many residues in %s !!!\n",AAdatafilename);
      break;
    }
    
    if ((strncmp(line,"END",3)) == 0) break;
  }

  fclose(AAdatafile);
}

AADatatable *LURAAData(char nameo[10],int create,AAdata dataofthisAminoAcid) {
  int h,i;

  h = hash(nameo);
  AADatatable *thiselement;

  for (thiselement = table[h]; thiselement != NULL; thiselement = thiselement->next) {
    if ((i=strcmp(nameo, thiselement->nameofthisresidue)) == 0) {
      return thiselement;
    }
  }
  if (create == 0){
    thiselement = (AADatatable *)emalloc(sizeof(AADatatable));
    for (i=0;nameo[i]!='\0';++i)
      thiselement->nameofthisresidue[i] = nameo[i];
    thiselement->nameofthisresidue[i]='\0';
    thiselement->AAdataontable = dataofthisAminoAcid;
    thiselement->next = table[h];
    table[h] = thiselement;
  }
  return thiselement;
}

unsigned int hash(char * str) {
  unsigned int h;
  unsigned char *p;

  h = 0;
  for (p = (unsigned char *)str; *p != '\0'; p++)
    h = MULTIPLIER * h + *p;
  h = h%NHASH;

  return h;
}
