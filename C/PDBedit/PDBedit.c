
#include <stdio.h>
#include <string.h>

#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "PDBedit.h"
#include "EF.h"
#include "PT.h"

int replace_ATOMNAME(PDBF PDB, char ATOMNAME1[5], char ATOMNAME2[5], int num){
  int i,len1,len;
  char ATOMNAMEdummy[4];

  len1=strlen(ATOMNAME1);
  len=strlen(ATOMNAME2);
  for (i=0;i<len1;++i) ATOMNAMEdummy[i]=PDB.PDBa[num].name[i+2];

  if ( strncmp(PDB.PDBa[num].name,ATOMNAME1,len1)==0) {
    for (i=0;i<len;++i) PDB.PDBa[num].name[2+i]=ATOMNAME2[i];
    for (i=len;i+2<5;++i) PDB.PDBa[num].name[2+i]=' ';
  }

}

int replaceall_ATOMNAME(PDBF PDB, char ATOMNAME[5], int num){
  int i,len;

  if ((len=strlen(ATOMNAME))>3) {
    printf("error:too long atom name\n");
    exit(1);
  }

  for (i=0;i<len;++i) PDB.PDBa[num].name[2+i]=ATOMNAME[i];
  for (i=len;i+2<5;++i) PDB.PDBa[num].name[2+i]=' ';

}

int delete_ATOMNAME(PDBF PDB, char ATOMNAME[5], int num){
  int i,j,numatom;
  int len1;
  char ATOMNAMEdummy[3];

  numatom=PDB.numatom;

  if ((len1=strlen(ATOMNAME))>3) {
    printf("error:too long atom name\n");
    exit(1);
  }

  for (i=0;i<3;++i) ATOMNAMEdummy[i]=PDB.PDBa[num].name[i+2];

  if ( strncmp(ATOMNAMEdummy,ATOMNAME,len1)==0) {
    for (i=num+1;i<numatom;++i) {
      PDB.PDBa[i-1].serial=PDB.PDBa[i-1].serial;        /* 7-11(serial)*/  
      for (j=0;j<4;++j)
	PDB.PDBa[i-1].name[j]=PDB.PDBa[i].name[j];      /*13-16(name)*/
      PDB.PDBa[i-1].altLOC=PDB.PDBa[i].altLOC;         /*17(alternate location indicator)*/
      for (j=0;j<3;++j)
	PDB.PDBa[i-1].resname[j]=PDB.PDBa[i].resname[j];   /*18-20(residue name)*/
      PDB.PDBa[i-1].ChainID=PDB.PDBa[i].ChainID;        /*22(Chain ID)*/
      PDB.PDBa[i-1].resSeq=PDB.PDBa[i].resSeq;        /*23-26(res seq)*/
      PDB.PDBa[i-1].iCode=PDB.PDBa[i].iCode;          /*27(iCode),28,29,30*/
      for (j=0;j<3;++j)                                 /*31-54(x,y,z)*/
	PDB.PDBa[i-1].coord[j]=PDB.PDBa[i].coord[j];
      PDB.PDBa[i-1].occupancy=PDB.PDBa[i].occupancy;  /*55-60(occupancy)*/
      PDB.PDBa[i-1].tempfact=PDB.PDBa[i].tempfact; /*61-66(tempfact)*/
    }

    //    PDB.numatom-=1;
    return 1;
  }

  return 0;
}

int pickup_ATOMNAME(PDBF PDB, char ATOMNAME[5],int num){
  int i,j,numatom;
  int len1;
  char ATOMNAMEdummy[3];

  numatom=PDB.numatom;

  if ((len1=strlen(ATOMNAME))>3) {
    printf("error:too long atom name\n");
    exit(1);
  }

  for (i=0;i<3;++i) ATOMNAMEdummy[i]=PDB.PDBa[num].name[i+2];

  if ( strncmp(ATOMNAMEdummy,ATOMNAME,len1)!=0) {
    for (i=num+1;i<numatom;++i) {
      PDB.PDBa[i-1].serial=PDB.PDBa[i-1].serial;        /* 7-11(serial)*/  
      for (j=0;j<4;++j)
	PDB.PDBa[i-1].name[j]=PDB.PDBa[i].name[j];      /*13-16(name)*/
      PDB.PDBa[i-1].altLOC=PDB.PDBa[i].altLOC;         /*17(alternate location indicator)*/
      for (j=0;j<3;++j)
	PDB.PDBa[i-1].resname[j]=PDB.PDBa[i].resname[j];   /*18-20(residue name)*/
      PDB.PDBa[i-1].ChainID=PDB.PDBa[i].ChainID;        /*22(Chain ID)*/
      PDB.PDBa[i-1].resSeq=PDB.PDBa[i].resSeq;        /*23-26(res seq)*/
      PDB.PDBa[i-1].iCode=PDB.PDBa[i].iCode;          /*27(iCode),28,29,30*/
      for (j=0;j<3;++j)                                 /*31-54(x,y,z)*/
	PDB.PDBa[i-1].coord[j]=PDB.PDBa[i].coord[j];
      PDB.PDBa[i-1].occupancy=PDB.PDBa[i].occupancy;  /*55-60(occupancy)*/
      PDB.PDBa[i-1].tempfact=PDB.PDBa[i].tempfact; /*61-66(tempfact)*/
    }
    //    PDB.numatom-=1;

    return 1;
  }

  return 0;
}

int put_occupancy(PDBF PDB,double occupancy,int num){
  int i,j;

  PDB.PDBa[num].occupancy=occupancy;
}

int show_line(PDBF PDB, int num){
  int i,j,numatom;

  printf("%s  ","ATOM");                   /* 1-6(ATOM  )*/
  printf("%5d",PDB.PDBa[num].serial);        /* 7-11(serial)*/  
  printf(" ");                             /* 12 */  
  for (j=0;j<4;++j)
    printf("%c",PDB.PDBa[num].name[j]);      /*13-16(name)*/
  printf("%c",PDB.PDBa[num].altLOC);         /*17(alternate location indicator)*/
  for (j=0;j<3;++j)
    printf("%c",PDB.PDBa[num].resname[j]);   /*18-20(residue name)*/
  printf(" ");                             /* 21 */  
  printf("%c",PDB.PDBa[num].ChainID);        /*22(Chain ID)*/
  printf("%4d",PDB.PDBa[num].resSeq);        /*23-26(res seq)*/
  printf("%c   ",PDB.PDBa[num].iCode);          /*27(iCode),28,29,30*/
  for (j=0;j<3;++j)                                 /*31-54(x,y,z)*/
    printf("%8.3lf",PDB.PDBa[num].coord[j]);
  printf("%6.2lf",PDB.PDBa[num].occupancy);  /*55-60(occupancy)*/
  printf("%6.2lf\n",PDB.PDBa[num].tempfact); /*61-66(tempfact)*/

}

int show_line_by_ATOMNAME(PDBF PDB,char ATOMNAME[5], int num){
  int i,j;

  if ( strncmp(PDB.PDBa[num].name,ATOMNAME,4)==0) {
    printf("%s  ","ATOM");                   /* 1-6(ATOM  )*/
    printf("%5d",PDB.PDBa[num].serial);        /* 7-11(serial)*/  
    printf(" ");                             /* 12 */  
    for (j=0;j<4;++j)
    printf ("%c",PDB.PDBa[num].name[j]);      /*13-16(name)*/
    printf("%c",PDB.PDBa[num].altLOC);         /*17(alternate location indicator)*/
    for (j=0;j<3;++j)
      printf("%c",PDB.PDBa[num].resname[j]);   /*18-20(residue name)*/
    printf(" ");                             /* 21 */  
    printf("%c",PDB.PDBa[num].ChainID);        /*22(Chain ID)*/
    printf("%4d",PDB.PDBa[num].resSeq);        /*23-26(res seq)*/
    printf("%c   ",PDB.PDBa[num].iCode);          /*27(iCode),28,29,30*/
    for (j=0;j<3;++j)                                 /*31-54(x,y,z)*/
      printf("%8.3lf",PDB.PDBa[num].coord[j]);
    printf("%6.2lf",PDB.PDBa[num].occupancy);  /*55-60(occupancy)*/
    printf("%6.2lf\n",PDB.PDBa[num].tempfact); /*61-66(tempfact)*/
  }
}

int show_line_by_RES(PDBF PDB,char *RESNAME, int num){
  int i,j;

  if ( strncmp(PDB.PDBa[num].resname,RESNAME,3)==0) {
    printf("%s  ","ATOM");                   /* 1-6(ATOM  )*/
    printf("%5d",PDB.PDBa[num].serial);        /* 7-11(serial)*/  
    printf(" ");                             /* 12 */  
    for (j=0;j<4;++j)
    printf("%c",PDB.PDBa[num].name[j]);      /*13-16(name)*/
    printf("%c",PDB.PDBa[num].altLOC);         /*17(alternate location indicator)*/
    for (j=0;j<3;++j)
      printf("%c",PDB.PDBa[num].resname[j]);   /*18-20(residue name)*/
    printf(" ");                             /* 21 */  
    printf("%c",PDB.PDBa[num].ChainID);        /*22(Chain ID)*/
    printf("%4d",PDB.PDBa[num].resSeq);        /*23-26(res seq)*/
    printf("%c   ",PDB.PDBa[num].iCode);          /*27(iCode),28,29,30*/
    for (j=0;j<3;++j)                                 /*31-54(x,y,z)*/
      printf("%8.3lf",PDB.PDBa[num].coord[j]);
    printf("%6.2lf",PDB.PDBa[num].occupancy);  /*55-60(occupancy)*/
    printf("%6.2lf\n",PDB.PDBa[num].tempfact); /*61-66(tempfact)*/
  }

}
