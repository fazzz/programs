#include <stdio.h>
#include <stdlib.h>

#include "PDB.h"

int writePDB(char *pdbfilenamebase, int flag) {  
  int i,j,k;
  int numatom;
  char pdbfilename[100];
  FILE *pdbfile;

  sprintf(pdbfilename,"%s.pdb",pdbfilenamebase);
  if((pdbfile=fopen(pdbfilename,flag))==NULL){
    printf("error: can not open %s \n",pdbfilename);
    exit(1);
  }

  numatom=PD.na;  
  for (i=0;i<numatom;++i) {
    for (j=0;j<5;++j)
      fprintf(pdbfile,"%c",PD.recordname[i][j]);/* 1-6(ATOM  )*/
    fprintf(pdbfile,"   %4d",PD.serial[i]);      /* 7-11(serial)*/  
    for (j=0;j<4;++j)
      fprintf(pdbfile,"%c",PD.name[i][j]);      /*13-16(name)*/
    fprintf(pdbfile,"%c ",PD.altLOC[i]);         /*17(alternate location indicator)*/
    for (j=0;j<3;++j)
      fprintf(pdbfile,"%c",PD.resname[i][j]);   /*18-20(residue name)*/
    fprintf(pdbfile,"%c  


",PD.ChainID[i]);      /*22(Chain ID)*/
    fprintf(pdbfile,"%4d",PD.resSeq[i]);        /*23-26(res seq)*/
    fprintf(pdbfile,"%c   ",PD.iCode[i]);      /*27(iCode)*/
    for (j=0;j<3;++j)                           /*31-54(x,y,z)*/
      fprintf(pdbfile,"%8.3lf",PD.coord[i][j]);
    fprintf(pdbfile,"%6.2lf",PD.occupancy[i]);  /*55-60(occupancy)*/
    fprintf(pdbfile,"%6.2lf\n",PD.tempfact[i]); /*61-66(tempfact)*/
  }
  return 1;
}
