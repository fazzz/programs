
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "PDB.h"
#include "EF.h"
#include "PT.h"

int writPDB(FILE *pdbfile,PDBF PDB) {  
  int i,j,k;
  int numatom;

  numatom=PDB.numatom;
  for (i=0;i<numatom;++i) {
    fprintf(pdbfile,"%s  ","ATOM");                   /* 1-6(ATOM  )*/
    fprintf(pdbfile,"%5d",PDB.PDBa[i].serial);        /* 7-11(serial)*/  
    fprintf(pdbfile," ");                             /* 12 */  
    for (j=0;j<4;++j)
      fprintf(pdbfile,"%c",PDB.PDBa[i].name[j]);      /*13-16(name)*/
    fprintf(pdbfile,"%c",PDB.PDBa[i].altLOC);         /*17(alternate location indicator)*/
    for (j=0;j<3;++j)
      fprintf(pdbfile,"%c",PDB.PDBa[i].resname[j]);   /*18-20(residue name)*/
    fprintf(pdbfile," ");                             /* 21 */  
    fprintf(pdbfile,"%c",PDB.PDBa[i].ChainID);        /*22(Chain ID)*/
    fprintf(pdbfile,"%4d",PDB.PDBa[i].resSeq);        /*23-26(res seq)*/
    fprintf(pdbfile,"%c   ",PDB.PDBa[i].iCode);          /*27(iCode),28,29,30*/
    for (j=0;j<3;++j)                                 /*31-54(x,y,z)*/
      fprintf(pdbfile,"%8.3lf",PDB.PDBa[i].coord[j]);
    fprintf(pdbfile,"%6.2lf",PDB.PDBa[i].occupancy);  /*55-60(occupancy)*/
    fprintf(pdbfile,"%6.2lf\n",PDB.PDBa[i].tempfact); /*61-66(tempfact)*/
  }
  return 1;
}

int writPDB_wopt(FILE *pdbfile,PDBF PDB, int MODE) {  
  int i,j,k;
  int numatom;
  int flag=ON;

  numatom=PDB.numatom;
  for (i=0;i<numatom;++i) {
    if (MODE == CA && strncmp(PDB.PDBa[i].name," CA",3)==0) flag=ON;
    else flag=OFF;
    if (MODE == HV && strncmp(PDB.PDBa[i].name," H",2)==0) flag=OFF;
    else flag=ON;
    if (MODE==AA) flag=ON;
    if (flag == ON) {
      fprintf(pdbfile,"%s  ","ATOM");                   /* 1-6(ATOM  )*/
      fprintf(pdbfile,"%5d",PDB.PDBa[i].serial);        /* 7-11(serial)*/  
      fprintf(pdbfile," ");                             /* 12 */  
      for (j=0;j<4;++j)
	fprintf(pdbfile,"%c",PDB.PDBa[i].name[j]);      /*13-16(name)*/
      fprintf(pdbfile,"%c",PDB.PDBa[i].altLOC);         /*17(alternate location indicator)*/
      for (j=0;j<3;++j)
	fprintf(pdbfile,"%c",PDB.PDBa[i].resname[j]);   /*18-20(residue name)*/
      fprintf(pdbfile," ");                             /* 21 */  
      fprintf(pdbfile,"%c",PDB.PDBa[i].ChainID);        /*22(Chain ID)*/
      fprintf(pdbfile,"%4d",PDB.PDBa[i].resSeq);        /*23-26(res seq)*/
      fprintf(pdbfile,"%c   ",PDB.PDBa[i].iCode);          /*27(iCode),28,29,30*/
      for (j=0;j<3;++j)                                 /*31-54(x,y,z)*/
	fprintf(pdbfile,"%8.3lf",PDB.PDBa[i].coord[j]);
      fprintf(pdbfile,"%6.2lf",PDB.PDBa[i].occupancy);  /*55-60(occupancy)*/
      fprintf(pdbfile,"%6.2lf\n",PDB.PDBa[i].tempfact); /*61-66(tempfact)*/
    }
  }
  return 1;
}

int readPDBdatafromParmtop(PDBF PDB) {
  int numatom,numres=-1,i;

  PDB.numatom=AP.NATOM;
  numatom=PDB.numatom;

  for (i=0;i<numatom;++i) {
    if (i==AP.IPRES[numres+1]-1) 
      ++numres;
    PDB.PDBa[i].HETEROflag=0;
    PDB.PDBa[i].serial=i+1;
    PDB.PDBa[i].name[0]=' '/*AP.IGRAPH[i][0]*/;
    PDB.PDBa[i].name[1]=AP.IGRAPH[i][0];
    PDB.PDBa[i].name[2]=AP.IGRAPH[i][1]/*' '*//*AP.IGRAPH[i][2]*/;
    PDB.PDBa[i].name[3]=AP.IGRAPH[i][2]/*' '*//*AP.IGRAPH[i][3]*/;
    PDB.PDBa[i].altLOC=' ';
    PDB.PDBa[i].resname[0]=/*' '*/AP.LABERES[numres][0];
    PDB.PDBa[i].resname[1]=/*' '*/AP.LABERES[numres][1];
    PDB.PDBa[i].resname[2]=/*' '*/AP.LABERES[numres][2];
    PDB.PDBa[i].ChainID=' ';
    PDB.PDBa[i].resSeq=/*0*/numres+1; // 11-2011
    PDB.PDBa[i].iCode=' ';
    PDB.PDBa[i].occupancy=0.0;
    PDB.PDBa[i].tempfact=0.0;
  }
}

int readPDBatomnum(FILE *pdbfile,int *numatom) {
  int c,n,na;
  int ATOMDATAFLAG=0;
  char recordname[6];

  na=0;
  n=-1;
  while ((c=getc(pdbfile))!=-1){
    if (c=='\n')
      n=-1;
    else
      ++n;
    if (n >= 0 && n < 6 )
      recordname[n]=c;
    if (n==6)
      if (strncmp(recordname,"ATOM",4)==0)
    	ATOMDATAFLAG=1;
      else if (strncmp(recordname,"HETATM",6)==0)
    	ATOMDATAFLAG=1;
    if (ATOMDATAFLAG==1 && n==6)
      ++na;
  }
  fclose(pdbfile);
  *numatom=na;
  return *numatom;
}

int readPDB(FILE *pdbfile,PDBF PDB,int numatom) {
  int i,j,c,n,na,d;
  int ATOMDATAFLAG=0;
  char recordname[6];
  double f,pmflag;

  na=-1;
  n=-1;
  
  while ((c=getc(pdbfile))!=-1){
    if (c=='\n')
      n=-1;
    else
      ++n;
    if (n >= 0 && n < 6 )
      recordname[n]=c;
    if (n==6)
      if (strncmp(recordname,"ATOM",4)==0)
    	ATOMDATAFLAG=1;
      else if (strncmp(recordname,"HETATM",6)==0)
    	ATOMDATAFLAG=1;
    if (ATOMDATAFLAG==1 && n==6)
      ++na;

    if (ATOMDATAFLAG==1){
      if (6 <= n && n <11) {
	if (c==' ')
	  d=0;
	else if (isdigit(c))
	  d=(c-'0')*pow(10,(11-n-1));
	PDB.PDBa[na].serial+=d;
      }
      else if (11<= n && n <16)
	PDB.PDBa[na].name[n-11]=c;
      else if (n == 16)
	PDB.PDBa[na].altLOC=c;
      else if (17<= n && n <20)
	PDB.PDBa[na].resname[n-17]=c;
      else if (n == 21)
	PDB.PDBa[na].ChainID=c;
      else if (22<= n && n <26) {
	if (c==' ')
	  d=0;
	else if (isdigit(c))
	  d=(c-'0')*pow(10,(26-n-1));
	PDB.PDBa[na].resSeq+=d;
      }
      else if (n==26)
	PDB.PDBa[na].iCode=c;
      else if (n==30)
	pmflag=1.0;
      else if (30<= n && n <34) {
	if (c==' ')
	  f=0;
	else if (c=='-')
	  pmflag=-1.0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(34-n-1));
	PDB.PDBa[na].coord[0]+=f;
      }
      else if (35<= n && n <38) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(35-n-1));
	PDB.PDBa[na].coord[0]+=f;
      }
      else if (n==37)
	  PDB.PDBa[na].coord[0]=pmflag*PDB.PDBa[na].coord[0];
      else if (n==38)
	pmflag=1.0;
      else if (38<= n && n <42) {
	if (c==' ')
	  f=0;
	else if (c=='-')
	  pmflag=-1.0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(42-n-1));
	PDB.PDBa[na].coord[1]+=f;
      }
      else if (43<= n && n <46) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(43-n-1));
	PDB.PDBa[na].coord[1]+=f;
      }
      else if (n==46)
	PDB.PDBa[na].coord[1]=pmflag*PDB.PDBa[na].coord[1];
      else if (n==46)
	pmflag=1.0;
      else if (46<= n && n <50) {
	if (c==' ')
	  f=0;
	else if (c=='-')
	  pmflag=-1.0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(50-n-1));
	PDB.PDBa[na].coord[2]+=f;
      }
      else if (51<= n && n <54) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(51-n-1));
	PDB.PDBa[na].coord[2]+=f;
      }
      else if (n==54) {
	PDB.PDBa[na].coord[2]=pmflag*PDB.PDBa[na].coord[2];
	PDB.PDBa[na].occupancy=0.0;
      }
      //      PDB.PDBa[na].occupancy=0.0;
      else if (54<= n && n <57) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(57-n-1));
	PDB.PDBa[na].occupancy+=f;
      }
      else if (58<= n && n <60) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(58-n-1));
	PDB.PDBa[na].occupancy+=f;
	PDB.PDBa[na].tempfact=0.0;
      }
      //      PDB.PDBa[na].tempfact=0.0;
      else if (60<= n && n <63) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(63-n-1));
	PDB.PDBa[na].tempfact+=f;
      }
      else if (64<= n && n <=66) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(64-n-1));
	PDB.PDBa[na].tempfact+=f;
      }
    }
  }
  //  PDB.numatom=numatom;
  fclose(pdbfile);
}

int readPDB2(FILE *pdbfile,PDBF PDB,int *numatom) {
  int i,j,c,n,na,d;
  int ATOMDATAFLAG=0;
  char recordname[6];
  double f,pmflag;

  na=-1;
  n=-1;
  
  while ((c=getc(pdbfile))!=-1){
    if (c=='\n')
      n=-1;
    else
      ++n;
    if (n >= 0 && n < 6 )
      recordname[n]=c;
    if (n==6)
      if (strncmp(recordname,"ATOM",4)==0)
    	ATOMDATAFLAG=1;
      else if (strncmp(recordname,"HETATM",6)==0)
    	ATOMDATAFLAG=1;
    if (ATOMDATAFLAG==1 && n==6)
      ++na;

    if (ATOMDATAFLAG==1){
      if (6 <= n && n <11) {
	if (c==' ')
	  d=0;
	else if (isdigit(c))
	  d=(c-'0')*pow(10,(11-n-1));
	PDB.PDBa[na].serial+=d;
      }
      if (11<= n && n <16)
	PDB.PDBa[na].name[n-11]=c;
      if (n == 16)
	PDB.PDBa[na].altLOC=c;
      if (17<= n && n <20)
	PDB.PDBa[na].resname[n-17]=c;
      if (n == 21)
	PDB.PDBa[na].ChainID=c;
      if (22<= n && n <26) {
	if (c==' ')
	  d=0;
	else if (isdigit(c))
	  d=(c-'0')*pow(10,(26-n-1));
	PDB.PDBa[na].resSeq+=d;
      }
      if (n==26)
	PDB.PDBa[na].iCode=c;
      if (n==30)
	pmflag=1.0;
      if (30<= n && n <34) {
	if (c==' ')
	  f=0;
	else if (c=='-')
	  pmflag=-1.0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(34-n-1));
	PDB.PDBa[na].coord[0]+=f;
      }
      if (35<= n && n <38) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(35-n-1));
	PDB.PDBa[na].coord[0]+=f;
      }
      if (n==37)
	  PDB.PDBa[na].coord[0]=pmflag*PDB.PDBa[na].coord[0];
      if (n==38)
	pmflag=1.0;
      if (38<= n && n <42) {
	if (c==' ')
	  f=0;
	else if (c=='-')
	  pmflag=-1.0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(42-n-1));
	PDB.PDBa[na].coord[1]+=f;
      }
      if (43<= n && n <46) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(43-n-1));
	PDB.PDBa[na].coord[1]+=f;
      }
      if (n==46)
	PDB.PDBa[na].coord[1]=pmflag*PDB.PDBa[na].coord[1];
      if (n==46)
	pmflag=1.0;
      if (46<= n && n <50) {
	if (c==' ')
	  f=0;
	else if (c=='-')
	  pmflag=-1.0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(50-n-1));
	PDB.PDBa[na].coord[2]+=f;
      }
      if (51<= n && n <54) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(51-n-1));
	PDB.PDBa[na].coord[2]+=f;
      }
      if (n==54)
	PDB.PDBa[na].coord[2]=pmflag*PDB.PDBa[na].coord[2];
      PDB.PDBa[na].occupancy=0.0;
      if (54<= n && n <57) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(57-n-1));
	PDB.PDBa[na].occupancy+=f;
      }
      if (58<= n && n <60) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(58-n-1));
	PDB.PDBa[na].occupancy+=f;
      }
      PDB.PDBa[na].tempfact=0.0;
      if (60<= n && n <63) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(63-n-1));
	PDB.PDBa[na].tempfact+=f;
      }
      if (64<= n && n <=66) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(64-n-1));
	PDB.PDBa[na].tempfact+=f;
      }
    }
  }
  *numatom=na;
  PDB.numatom=*numatom;
  fclose(pdbfile);
}

int copyPDBform(PDBF PDB1,PDBF PDB2) {
  int i,j,numatom;

  PDB2.numatom=PDB1.numatom;
  numatom=PDB2.numatom;

  for (i=0;i<numatom;++i) {
    PDB2.PDBa[i].serial    =PDB1.PDBa[i].serial;    
    PDB2.PDBa[i].name[0]   =PDB1.PDBa[i].name[0];   
    PDB2.PDBa[i].name[1]   =PDB1.PDBa[i].name[1];   
    PDB2.PDBa[i].name[2]   =PDB1.PDBa[i].name[2];   
    PDB2.PDBa[i].name[3]   =PDB1.PDBa[i].name[3];   
    PDB2.PDBa[i].altLOC    =PDB1.PDBa[i].altLOC;    
    PDB2.PDBa[i].resname[0]=PDB1.PDBa[i].resname[0];
    PDB2.PDBa[i].resname[1]=PDB1.PDBa[i].resname[1];
    PDB2.PDBa[i].resname[2]=PDB1.PDBa[i].resname[2];
    PDB2.PDBa[i].resname[3]=PDB1.PDBa[i].resname[3];
    PDB2.PDBa[i].ChainID   =PDB1.PDBa[i].ChainID;   
    PDB2.PDBa[i].resSeq    =PDB1.PDBa[i].resSeq;    
    PDB2.PDBa[i].iCode     =PDB1.PDBa[i].iCode;     
    for (j=0;j<3;++j)
      PDB2.PDBa[i].coord[j]   =PDB1.PDBa[i].coord[j];
    PDB2.PDBa[i].occupancy =PDB1.PDBa[i].occupancy; 
    PDB2.PDBa[i].tempfact  =PDB1.PDBa[i].tempfact;  
  }
}

int addCAP(PDBF PDB1,PDBF PDB2) {
  int i,j,numatom;
  double cN[3],cCA[3],vCAN[3];
  double cC[3],cO[3],vCO[3];

  double lCC=1.52200000;

  numatom=PDB2.numatom;

  PDB2.PDBa[0].serial    =PDB1.PDBa[i-1].serial;    

  for (i=0;i<PDB1.numatom;++i) {
    if (strncmp(PDB1.PDBa[i].name,"  N ",4)==0) {
      for (j=0;j<3;++j)	cN[j]=PDB1.PDBa[i].coord[j];
    }
    if (strncmp(PDB1.PDBa[i].name,"  CA",4)==0) {
      for (j=0;j<3;++j)	cCA[j]=PDB1.PDBa[i].coord[j];
      break;
    }
  }

  PDB2.PDBa[0].serial=0;
  PDB2.PDBa[0].name[0]   =' ';   
  PDB2.PDBa[0].name[1]   =' ';   
  PDB2.PDBa[0].name[2]   ='C';   
  PDB2.PDBa[0].name[3]   =' ';
  PDB2.PDBa[0].altLOC     =' ';
  PDB2.PDBa[0].resname[0] ='A';
  PDB2.PDBa[0].resname[1] ='C';
  PDB2.PDBa[0].resname[2] ='E';
  PDB2.PDBa[0].ChainID=' ';        
  PDB2.PDBa[0].resSeq=0;        
  PDB2.PDBa[0].iCode=' ';       
  for (j=0;j<3;++j)                                 
    PDB2.PDBa[0].coord[j]=cN[j]+lCC*(cN[j]-cCA[j]);
  PDB2.PDBa[0].occupancy=1.0;  
  PDB2.PDBa[0].tempfact=0.0; 

  for (i=1;i<numatom-1;++i) {
    PDB2.PDBa[i].serial    =PDB1.PDBa[i-1].serial;    
    PDB2.PDBa[i].name[0]   =PDB1.PDBa[i-1].name[0];   
    PDB2.PDBa[i].name[1]   =PDB1.PDBa[i-1].name[1];   
    PDB2.PDBa[i].name[2]   =PDB1.PDBa[i-1].name[2];   
    PDB2.PDBa[i].name[3]   =PDB1.PDBa[i-1].name[3];   
    PDB2.PDBa[i].altLOC    =PDB1.PDBa[i-1].altLOC;    
    PDB2.PDBa[i].resname[0]=PDB1.PDBa[i-1].resname[0];
    PDB2.PDBa[i].resname[1]=PDB1.PDBa[i-1].resname[1];
    PDB2.PDBa[i].resname[2]=PDB1.PDBa[i-1].resname[2];
    PDB2.PDBa[i].resname[3]=PDB1.PDBa[i-1].resname[3];
    PDB2.PDBa[i].ChainID   =PDB1.PDBa[i-1].ChainID;   
    PDB2.PDBa[i].resSeq    =PDB1.PDBa[i-1].resSeq;    
    PDB2.PDBa[i].iCode     =PDB1.PDBa[i-1].iCode;     
    for (j=0;j<3;++j)
      PDB2.PDBa[i].coord[j]   =PDB1.PDBa[i-1].coord[j];
    PDB2.PDBa[i].occupancy =PDB1.PDBa[i-1].occupancy; 
    PDB2.PDBa[i].tempfact  =PDB1.PDBa[i-1].tempfact;  
  }

  for (i=PDB1.numatom-1;i>=0;--i) {
    if (strncmp(PDB1.PDBa[i].name,"  O ",4)==0) {
      for (j=0;j<3;++j)	cO[j]=PDB1.PDBa[i].coord[j];
    }
    if (strncmp(PDB1.PDBa[i].name,"  C ",4)==0) {
      for (j=0;j<3;++j)	cC[j]=PDB1.PDBa[i].coord[j];
      break;
    }
  }

  PDB2.PDBa[PDB1.numatom].serial=PDB1.numatom+1;
  PDB2.PDBa[PDB1.numatom].name[0]   =' ';   
  PDB2.PDBa[PDB1.numatom].name[1]   =' ';   
  PDB2.PDBa[PDB1.numatom].name[2]   ='N';   
  PDB2.PDBa[PDB1.numatom].name[3]   =' ';
  PDB2.PDBa[PDB1.numatom].altLOC     =' ';
  PDB2.PDBa[PDB1.numatom].resname[0] ='N';
  PDB2.PDBa[PDB1.numatom].resname[1] ='M';
  PDB2.PDBa[PDB1.numatom].resname[2] ='E';
  PDB2.PDBa[PDB1.numatom].ChainID=' ';        
  PDB2.PDBa[PDB1.numatom].resSeq=PDB2.PDBa[PDB1.numatom-1].resSeq+1;
  PDB2.PDBa[PDB1.numatom].iCode=' ';       
  for (j=0;j<3;++j)                                 
    PDB2.PDBa[PDB1.numatom].coord[j]=cO[j]+lCC*(cO[j]-cC[j]);
  PDB2.PDBa[PDB1.numatom].occupancy=1.0;  
  PDB2.PDBa[PDB1.numatom].tempfact=0.0; 

}
