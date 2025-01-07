
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include "PDBL.h"
#include "EF.h"
#include "PTL.h"

int writPDBL(FILE *pdbfile,PDBLF PDBL) {  
  int i,j,k;
  int numatom;

  numatom=PDBL.numatom;
  for (i=0;i<numatom;++i) {
    fprintf(pdbfile,"%s  ","ATOM");                   /* 1-6(ATOM  )*/
    fprintf(pdbfile,"%5d",PDBL.PDBLa[i].serial);        /* 7-11(serial)*/  
    fprintf(pdbfile," ");                             /* 12 */  
    for (j=0;j<4;++j)
      fprintf(pdbfile,"%c",PDBL.PDBLa[i].name[j]);      /*13-16(name)*/
    fprintf(pdbfile,"%c",PDBL.PDBLa[i].altLOC);         /*17(alternate location indicator)*/
    for (j=0;j<3;++j)
      fprintf(pdbfile,"%c",PDBL.PDBLa[i].resname[j]);   /*18-20(residue name)*/
    fprintf(pdbfile," ");                             /* 21 */  
    fprintf(pdbfile,"%c",PDBL.PDBLa[i].ChainID);        /*22(Chain ID)*/
    fprintf(pdbfile,"%4d",PDBL.PDBLa[i].resSeq);        /*23-26(res seq)*/
    fprintf(pdbfile,"%c   ",PDBL.PDBLa[i].iCode);          /*27(iCode),28,29,30*/
    for (j=0;j<3;++j)                                 /*31-54(x,y,z)*/
      fprintf(pdbfile,"%8.3lf",PDBL.PDBLa[i].coord[j]);
    fprintf(pdbfile,"%6.2lf",PDBL.PDBLa[i].occupancy);  /*55-60(occupancy)*/
    fprintf(pdbfile,"%6.2lf\n",PDBL.PDBLa[i].tempfact); /*61-66(tempfact)*/
  }
  return 1;
}

int writPDBL_wopt(FILE *pdbfile,PDBLF PDBL, int MODE) {  
  int i,j,k;
  int numatom;
  int flag=ON;

  numatom=PDBL.numatom;
  for (i=0;i<numatom;++i) {
    if (MODE == CA ) {
      if (strncmp(PDBL.PDBLa[i].name," CA",3)==0) flag=ON;
      else flag=OFF;
    }
    else if (MODE == HV ) {
      if ( strncmp(PDBL.PDBLa[i].name," H",2)==0 ) flag=OFF;
      else flag=ON;
    }
    else if (MODE==AA) flag=ON;

    if (flag == ON) {
      fprintf(pdbfile,"%s  ","ATOM");                   /* 1-6(ATOM  )*/
      fprintf(pdbfile,"%5d",PDBL.PDBLa[i].serial);        /* 7-11(serial)*/  
      fprintf(pdbfile," ");                             /* 12 */  
      for (j=0;j<4;++j) {
	if (PDBL.PDBLa[i].name[j]!='\0')
	  fprintf(pdbfile,"%c",PDBL.PDBLa[i].name[j]);      /*13-16(name)*/
	else 
	  break;
      }
      for (;j<4;++j)
	fprintf(pdbfile," ",PDBL.PDBLa[i].name[j]);      /*13-16(name)*/
      fprintf(pdbfile,"%c",PDBL.PDBLa[i].altLOC);         /*17(alternate location indicator)*/
      for (j=0;j<3;++j) {
	if (PDBL.PDBLa[i].resname[j]!='\0')
	  fprintf(pdbfile,"%c",PDBL.PDBLa[i].resname[j]);   /*18-20(residue name)*/
	else 
	  break;
      }
      for (;j<4;++j)
	fprintf(pdbfile," ",PDBL.PDBLa[i].name[j]);      /*13-16(name)*/
      //      fprintf(pdbfile," ");                             /* 21 */  
      fprintf(pdbfile,"%c",PDBL.PDBLa[i].ChainID);        /*22(Chain ID)*/
      fprintf(pdbfile,"%4d",PDBL.PDBLa[i].resSeq);        /*23-26(res seq)*/
      fprintf(pdbfile,"%c   ",PDBL.PDBLa[i].iCode);          /*27(iCode),28,29,30*/
      for (j=0;j<3;++j)                                 /*31-54(x,y,z)*/
	fprintf(pdbfile,"%8.3lf",PDBL.PDBLa[i].coord[j]);
      fprintf(pdbfile,"%6.2lf",PDBL.PDBLa[i].occupancy);  /*55-60(occupancy)*/
      fprintf(pdbfile,"%6.2lf\n",PDBL.PDBLa[i].tempfact); /*61-66(tempfact)*/
    }
  }
  return 1;
}

int writPDBL_wopt_series(FILE *pdbfile,PDBLF PDBL, int MODE) {  
  int i,j,k;
  int numatom;
  int flag=ON;

  fprintf(pdbfile,"MODEL \n");
  numatom=PDBL.numatom;
  for (i=0;i<numatom;++i) {
    if (MODE == CA ) {
      if (strncmp(PDBL.PDBLa[i].name," CA",3)==0) flag=ON;
      else flag=OFF;
    }
    else if (MODE == HV ) {
      if ( strncmp(PDBL.PDBLa[i].name," H",2)==0 ) flag=OFF;
      else flag=ON;
    }
    else if (MODE==AA) flag=ON;

    if (flag == ON) {
      fprintf(pdbfile,"%s  ","ATOM");                   /* 1-6(ATOM  )*/
      fprintf(pdbfile,"%5d",PDBL.PDBLa[i].serial);        /* 7-11(serial)*/  
      fprintf(pdbfile," ");                             /* 12 */  
      for (j=0;j<4;++j) {
	if (PDBL.PDBLa[i].name[j]!='\0')
	  fprintf(pdbfile,"%c",PDBL.PDBLa[i].name[j]);      /*13-16(name)*/
	else 
	  break;
      }
      for (;j<4;++j)
	fprintf(pdbfile," ",PDBL.PDBLa[i].name[j]);      /*13-16(name)*/
      fprintf(pdbfile,"%c",PDBL.PDBLa[i].altLOC);         /*17(alternate location indicator)*/
      for (j=0;j<3;++j) {
	if (PDBL.PDBLa[i].resname[j]!='\0')
	  fprintf(pdbfile,"%c",PDBL.PDBLa[i].resname[j]);   /*18-20(residue name)*/
	else 
	  break;
      }
      for (;j<4;++j)
	fprintf(pdbfile," ",PDBL.PDBLa[i].name[j]);      /*13-16(name)*/
      //      fprintf(pdbfile," ");                             /* 21 */  
      fprintf(pdbfile,"%c",PDBL.PDBLa[i].ChainID);        /*22(Chain ID)*/
      fprintf(pdbfile,"%4d",PDBL.PDBLa[i].resSeq);        /*23-26(res seq)*/
      fprintf(pdbfile,"%c   ",PDBL.PDBLa[i].iCode);          /*27(iCode),28,29,30*/
      for (j=0;j<3;++j)                                 /*31-54(x,y,z)*/
	fprintf(pdbfile,"%8.3lf",PDBL.PDBLa[i].coord[j]);
      fprintf(pdbfile,"%6.2lf",PDBL.PDBLa[i].occupancy);  /*55-60(occupancy)*/
      fprintf(pdbfile,"%6.2lf\n",PDBL.PDBLa[i].tempfact); /*61-66(tempfact)*/
    }
  }
  fprintf(pdbfile,"ENDMOD\n");
  return 1;
}

int readPDBLdatafromParmtop(PDBLF PDBL) {
  int numatom,numres=-1,i;

  PDBL.numatom=AP.NATOM;
  numatom=PDBL.numatom;

  for (i=0;i<numatom;++i) {
    if (i==AP.IPRES[numres+1]-1) 
      ++numres;
    PDBL.PDBLa[i].HETEROflag=0;
    PDBL.PDBLa[i].serial=i+1;
    if (strncmp(AP.IGRAPH[i],"HH3",3)!=0) {
      PDBL.PDBLa[i].name[0]=' '/*AP.IGRAPH[i][0]*/;
      PDBL.PDBLa[i].name[1]=AP.IGRAPH[i][0];
      PDBL.PDBLa[i].name[2]=AP.IGRAPH[i][1]/*' '*//*AP.IGRAPH[i][2]*/;
      PDBL.PDBLa[i].name[3]=AP.IGRAPH[i][2]/*' '*//*AP.IGRAPH[i][3]*/;
    }
    else {
      PDBL.PDBLa[i].name[0]=AP.IGRAPH[i][0];
      PDBL.PDBLa[i].name[1]=AP.IGRAPH[i][1];
      PDBL.PDBLa[i].name[2]=AP.IGRAPH[i][2];
      PDBL.PDBLa[i].name[3]=AP.IGRAPH[i][3];
    }
    PDBL.PDBLa[i].altLOC=' ';
    PDBL.PDBLa[i].resname[0]=/*' '*/AP.LABERES[numres][0];
    PDBL.PDBLa[i].resname[1]=/*' '*/AP.LABERES[numres][1];
    PDBL.PDBLa[i].resname[2]=/*' '*/AP.LABERES[numres][2];
    PDBL.PDBLa[i].ChainID=' ';
    PDBL.PDBLa[i].resSeq=/*0*/numres+1; // 11-2011
    PDBL.PDBLa[i].iCode=' ';
    PDBL.PDBLa[i].occupancy=0.0;
    PDBL.PDBLa[i].tempfact=0.0;
  }
}

int readPDBLatomnum(FILE *pdbfile,int *numatom) {
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

int readPDBL(FILE *pdbfile,PDBLF PDBL,int numatom) {
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
	PDBL.PDBLa[na].serial+=d;
      }
      if (11<= n && n <16)
	PDBL.PDBLa[na].name[n-11]=c;
      if (n == 16)
	PDBL.PDBLa[na].altLOC=c;
      if (17<= n && n <20)
	PDBL.PDBLa[na].resname[n-17]=c;
      if (n == 21)
	PDBL.PDBLa[na].ChainID=c;
      if (22<= n && n <26) {
	if (c==' ')
	  d=0;
	else if (isdigit(c))
	  d=(c-'0')*pow(10,(26-n-1));
	PDBL.PDBLa[na].resSeq+=d;
      }
      if (n==26)
	PDBL.PDBLa[na].iCode=c;
      if (n==30)
	pmflag=1.0;
      if (30<= n && n <34) {
	if (c==' ')
	  f=0;
	else if (c=='-')
	  pmflag=-1.0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(34-n-1));
	PDBL.PDBLa[na].coord[0]+=f;
      }
      if (35<= n && n <38) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(35-n-1));
	PDBL.PDBLa[na].coord[0]+=f;
      }
      if (n==37)
	  PDBL.PDBLa[na].coord[0]=pmflag*PDBL.PDBLa[na].coord[0];
      if (n==38)
	pmflag=1.0;
      if (38<= n && n <42) {
	if (c==' ')
	  f=0;
	else if (c=='-')
	  pmflag=-1.0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(42-n-1));
	PDBL.PDBLa[na].coord[1]+=f;
      }
      if (43<= n && n <46) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(43-n-1));
	PDBL.PDBLa[na].coord[1]+=f;
      }
      if (n==46)
	PDBL.PDBLa[na].coord[1]=pmflag*PDBL.PDBLa[na].coord[1];
      if (n==46)
	pmflag=1.0;
      if (46<= n && n <50) {
	if (c==' ')
	  f=0;
	else if (c=='-')
	  pmflag=-1.0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(50-n-1));
	PDBL.PDBLa[na].coord[2]+=f;
      }
      if (51<= n && n <54) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(51-n-1));
	PDBL.PDBLa[na].coord[2]+=f;
      }
      if (n==54)
	PDBL.PDBLa[na].coord[2]=pmflag*PDBL.PDBLa[na].coord[2];
      PDBL.PDBLa[na].occupancy=0.0;
      if (54<= n && n <57) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(57-n-1));
	PDBL.PDBLa[na].occupancy+=f;
      }
      if (58<= n && n <60) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(58-n-1));
	PDBL.PDBLa[na].occupancy+=f;
      }
      PDBL.PDBLa[na].tempfact=0.0;
      if (60<= n && n <63) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(63-n-1));
	PDBL.PDBLa[na].tempfact+=f;
      }
      if (64<= n && n <=66) {
	if (c==' ')
	  f=0;
	else if (isdigit(c))
	  f=(double)(c-'0')*pow(10.0,(64-n-1));
	PDBL.PDBLa[na].tempfact+=f;
      }
    }
  }
  //  PDBL.numatom=numatom;
  fclose(pdbfile);
}

int copyPDBLform(PDBLF PDBL1,PDBLF PDBL2) {
  int i,j,numatom;

  PDBL2.numatom=PDBL1.numatom;
  numatom=PDBL2.numatom;

  for (i=0;i<numatom;++i) {
    PDBL2.PDBLa[i].serial    =PDBL1.PDBLa[i].serial;    
    PDBL2.PDBLa[i].name[0]   =PDBL1.PDBLa[i].name[0];   
    PDBL2.PDBLa[i].name[1]   =PDBL1.PDBLa[i].name[1];   
    PDBL2.PDBLa[i].name[2]   =PDBL1.PDBLa[i].name[2];   
    PDBL2.PDBLa[i].name[3]   =PDBL1.PDBLa[i].name[3];   
    PDBL2.PDBLa[i].altLOC    =PDBL1.PDBLa[i].altLOC;    
    PDBL2.PDBLa[i].resname[0]=PDBL1.PDBLa[i].resname[0];
    PDBL2.PDBLa[i].resname[1]=PDBL1.PDBLa[i].resname[1];
    PDBL2.PDBLa[i].resname[2]=PDBL1.PDBLa[i].resname[2];
    PDBL2.PDBLa[i].resname[3]=PDBL1.PDBLa[i].resname[3];
    PDBL2.PDBLa[i].ChainID   =PDBL1.PDBLa[i].ChainID;   
    PDBL2.PDBLa[i].resSeq    =PDBL1.PDBLa[i].resSeq;    
    PDBL2.PDBLa[i].iCode     =PDBL1.PDBLa[i].iCode;     
    for (j=0;j<3;++j)
      PDBL2.PDBLa[i].coord[j]   =PDBL1.PDBLa[i].coord[j];
    PDBL2.PDBLa[i].occupancy =PDBL1.PDBLa[i].occupancy; 
    PDBL2.PDBLa[i].tempfact  =PDBL1.PDBLa[i].tempfact;  
  }
}

int addCAPL(PDBLF PDBL1,PDBLF PDBL2) {
  int i,j,numatom;
  double cN[3],cCA[3],vCAN[3];
  double cC[3],cO[3],vCO[3];

  double lCC=1.52200000;

  numatom=PDBL2.numatom;

  PDBL2.PDBLa[0].serial    =PDBL1.PDBLa[i-1].serial;    

  for (i=0;i<PDBL1.numatom;++i) {
    if (strncmp(PDBL1.PDBLa[i].name,"  N ",4)==0) {
      for (j=0;j<3;++j)	cN[j]=PDBL1.PDBLa[i].coord[j];
    }
    if (strncmp(PDBL1.PDBLa[i].name,"  CA",4)==0) {
      for (j=0;j<3;++j)	cCA[j]=PDBL1.PDBLa[i].coord[j];
      break;
    }
  }

  PDBL2.PDBLa[0].serial=0;
  PDBL2.PDBLa[0].name[0]   =' ';   
  PDBL2.PDBLa[0].name[1]   =' ';   
  PDBL2.PDBLa[0].name[2]   ='C';   
  PDBL2.PDBLa[0].name[3]   =' ';
  PDBL2.PDBLa[0].altLOC     =' ';
  PDBL2.PDBLa[0].resname[0] ='A';
  PDBL2.PDBLa[0].resname[1] ='C';
  PDBL2.PDBLa[0].resname[2] ='E';
  PDBL2.PDBLa[0].ChainID=' ';        
  PDBL2.PDBLa[0].resSeq=0;        
  PDBL2.PDBLa[0].iCode=' ';       
  for (j=0;j<3;++j)                                 
    PDBL2.PDBLa[0].coord[j]=cN[j]+lCC*(cN[j]-cCA[j]);
  PDBL2.PDBLa[0].occupancy=1.0;  
  PDBL2.PDBLa[0].tempfact=0.0; 

  for (i=1;i<numatom-1;++i) {
    PDBL2.PDBLa[i].serial    =PDBL1.PDBLa[i-1].serial;    
    PDBL2.PDBLa[i].name[0]   =PDBL1.PDBLa[i-1].name[0];   
    PDBL2.PDBLa[i].name[1]   =PDBL1.PDBLa[i-1].name[1];   
    PDBL2.PDBLa[i].name[2]   =PDBL1.PDBLa[i-1].name[2];   
    PDBL2.PDBLa[i].name[3]   =PDBL1.PDBLa[i-1].name[3];   
    PDBL2.PDBLa[i].altLOC    =PDBL1.PDBLa[i-1].altLOC;    
    PDBL2.PDBLa[i].resname[0]=PDBL1.PDBLa[i-1].resname[0];
    PDBL2.PDBLa[i].resname[1]=PDBL1.PDBLa[i-1].resname[1];
    PDBL2.PDBLa[i].resname[2]=PDBL1.PDBLa[i-1].resname[2];
    PDBL2.PDBLa[i].resname[3]=PDBL1.PDBLa[i-1].resname[3];
    PDBL2.PDBLa[i].ChainID   =PDBL1.PDBLa[i-1].ChainID;   
    PDBL2.PDBLa[i].resSeq    =PDBL1.PDBLa[i-1].resSeq;    
    PDBL2.PDBLa[i].iCode     =PDBL1.PDBLa[i-1].iCode;     
    for (j=0;j<3;++j)
      PDBL2.PDBLa[i].coord[j]   =PDBL1.PDBLa[i-1].coord[j];
    PDBL2.PDBLa[i].occupancy =PDBL1.PDBLa[i-1].occupancy; 
    PDBL2.PDBLa[i].tempfact  =PDBL1.PDBLa[i-1].tempfact;  
  }

  for (i=PDBL1.numatom-1;i>=0;--i) {
    if (strncmp(PDBL1.PDBLa[i].name,"  O ",4)==0) {
      for (j=0;j<3;++j)	cO[j]=PDBL1.PDBLa[i].coord[j];
    }
    if (strncmp(PDBL1.PDBLa[i].name,"  C ",4)==0) {
      for (j=0;j<3;++j)	cC[j]=PDBL1.PDBLa[i].coord[j];
      break;
    }
  }

  PDBL2.PDBLa[PDBL1.numatom].serial=PDBL1.numatom+1;
  PDBL2.PDBLa[PDBL1.numatom].name[0]   =' ';   
  PDBL2.PDBLa[PDBL1.numatom].name[1]   =' ';   
  PDBL2.PDBLa[PDBL1.numatom].name[2]   ='N';   
  PDBL2.PDBLa[PDBL1.numatom].name[3]   =' ';
  PDBL2.PDBLa[PDBL1.numatom].altLOC     =' ';
  PDBL2.PDBLa[PDBL1.numatom].resname[0] ='N';
  PDBL2.PDBLa[PDBL1.numatom].resname[1] ='M';
  PDBL2.PDBLa[PDBL1.numatom].resname[2] ='E';
  PDBL2.PDBLa[PDBL1.numatom].ChainID=' ';        
  PDBL2.PDBLa[PDBL1.numatom].resSeq=PDBL2.PDBLa[PDBL1.numatom-1].resSeq+1;
  PDBL2.PDBLa[PDBL1.numatom].iCode=' ';       
  for (j=0;j<3;++j)                                 
    PDBL2.PDBLa[PDBL1.numatom].coord[j]=cO[j]+lCC*(cO[j]-cC[j]);
  PDBL2.PDBLa[PDBL1.numatom].occupancy=1.0;  
  PDBL2.PDBLa[PDBL1.numatom].tempfact=0.0; 

}
