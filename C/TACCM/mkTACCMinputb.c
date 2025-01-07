
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

#include "PTL.h"
#include "mkTACCMinputb.h"
#include "AminoAcid.h"

int mkTACCMinput_set_atomnumi_atomnuml(int atomnumj,int atomnumk,
				       char *atomnamej, char *atomnamek,
				       char *resnamej,  char *resnamek,
				       int *atomnumi,  int *atomnuml, int *mesg, int termflag) {
  int i,j,k;
  int numi=0,numl=0;
  int c;

  int flag=S1;
  
  char atomnamei[4],atomnamel[4];
  //  char treei[4],treel[4];

  *mesg=1;

  c=judge_dihedtype(atomnamej,atomnamek/*,treej,treek*/,resnamej,resnamek,
		    atomnamei,atomnamel/*,treei,treel*/,termflag);

  for (i=0;i<sizeof(atomnamei)/sizeof(char);++i) if (atomnamei[i]=='\0') break; numi=i;
  for (i=0;i<sizeof(atomnamei)/sizeof(char);++i) if (atomnamel[i]=='\0') break; numl=i;

  if (strncmp(atomnamei,"H",1)!=0) {
    for (i=0;i<AP.MBONA;++i) {
      if (abs(AP.BA[i][0])/3==atomnumj) {
	if (strncmp(AP.IGRAPH[abs(AP.BA[i][1])/3],atomnamei,numi)==0) {
	  if (abs(AP.BA[i][1])/3 != atomnumk && 
	      ( strncmp(AP.ITREE[abs(AP.BA[i][1])/3],"M",1)==0 || c > OMEGA || termflag != NOTERM ) ||
	      /*strncmp(AP.ITREE[abs(AP.BA[i][1])/3],"BLA",1)!=0*/ 
	      ( abs(AP.BA[i][1])/3 /*>=*/< AP.IPRES[1]-1 || abs(AP.BA[i][1])/3 >=/*<*/ AP.IPRES[AP.NRES-1]-1 ) ) {
	    *atomnumi=abs(AP.BA[i][1])/3;
	    flag=S2;
	    break;
	  }
	}
      }    
      else if (abs(AP.BA[i][1])/3==atomnumj) {
	if (strncmp(AP.IGRAPH[abs(AP.BA[i][0])/3],atomnamei,numi)==0) {
	  if (abs(AP.BA[i][0])/3 != atomnumk && 
	      ( strncmp(AP.ITREE[abs(AP.BA[i][0])/3],"M",1)==0 || c > OMEGA || termflag != NOTERM ) ||
	      /*strncmp(AP.ITREE[abs(AP.BA[i][1])/3],"BLA",1)!=0*/ 
	      ( abs(AP.BA[i][0])/3 /*>=*/< AP.IPRES[1]-1 || abs(AP.BA[i][1])/3 >=/*<*/ AP.IPRES[AP.NRES-1]-1 ) ) {
	    *atomnumi=abs(AP.BA[i][0])/3;
	    flag=S2;
	    break;
	  }
	}
      }
    }
  }
  else {
    for (i=0;i<AP.NBONH;++i) {
      if (abs(AP.BH[i][0])/3==atomnumj) {
	if (strncmp(AP.IGRAPH[abs(AP.BH[i][1])/3],atomnamei,numi)==0) {
	  if (abs(AP.BH[i][1])/3 != atomnumk  && 
	      ( strncmp(AP.ITREE[abs(AP.BH[i][1])/3],"M",1)==0 || c > OMEGA || termflag != NOTERM ) ||
	      /*strncmp(AP.ITREE[abs(AP.BA[i][1])/3],"BLA",1)!=0*/ 
	      ( abs(AP.BH[i][1])/3 /*>=*/< AP.IPRES[1]-1 || abs(AP.BH[i][1])/3 >=/*<*/ AP.IPRES[AP.NRES-1]-1 ) ) {
	    *atomnumi=abs(AP.BH[i][1])/3;
	    flag=S2;
	    break;
	  }
	}
      }    
      else if (abs(AP.BH[i][1])/3==atomnumj) {
	if (strncmp(AP.IGRAPH[abs(AP.BH[i][0])/3],atomnamei,numi)==0) {
	  if (abs(AP.BH[i][0])/3 != atomnumk  && 
	      ( strncmp(AP.ITREE[abs(AP.BH[i][0])/3],"M",1)==0 || c > OMEGA || termflag != NOTERM ) ||
	      /*strncmp(AP.ITREE[abs(AP.BA[i][1])/3],"BLA",1)!=0*/ 
	      ( abs(AP.BH[i][0])/3 /*>=*/< AP.IPRES[1]-1 || abs(AP.BH[i][0])/3 >=/*<*/ AP.IPRES[AP.NRES-1]-1 ) ) {
	    *atomnumi=abs(AP.BH[i][0])/3;
	    flag=S2;
	    break;
	  }
	}
      }
    }
  }

  if (strncmp(atomnamel,"H",1)!=0) {
    for (i=0;i<AP.MBONA;++i) {
      if (abs(AP.BA[i][0])/3==atomnumk) {
	if (strncmp(AP.IGRAPH[abs(AP.BA[i][1])/3],atomnamel,numl)==0) {
	  if (abs(AP.BA[i][1])/3 != atomnumj  && 
	      ( strncmp(AP.ITREE[abs(AP.BA[i][1])/3],"M",1)==0 || c > OMEGA || termflag != NOTERM ) ||
	      /*strncmp(AP.ITREE[abs(AP.BA[i][1])/3],"BLA",1)!=0*/ 
	      ( abs(AP.BA[i][1])/3 /*>=*/< AP.IPRES[1]-1 || abs(AP.BA[i][1])/3 /*<*/>= AP.IPRES[AP.NRES-1]-1 ) ) {
	    *atomnuml=abs(AP.BA[i][1])/3;
	    flag=S3;
	    break;
	  }
	}
      }    
      else if (abs(AP.BA[i][1])/3==atomnumk) {
	if (strncmp(AP.IGRAPH[abs(AP.BA[i][0])/3],atomnamel,numl)==0) {
	  if (abs(AP.BA[i][0])/3 != atomnumj  && 
	      ( strncmp(AP.ITREE[abs(AP.BA[i][0])/3],"M",1)==0 || c > OMEGA || termflag != NOTERM ) ||
	      /*strncmp(AP.ITREE[abs(AP.BA[i][1])/3],"BLA",1)!=0*/ 
	      ( abs(AP.BA[i][1])/3 /*>=*/< AP.IPRES[1]-1 || abs(AP.BA[i][1])/3 /*<*/>= AP.IPRES[AP.NRES-1]-1 ) ) {
	    *atomnuml=abs(AP.BA[i][0])/3;
	    flag=S3;
	    break;
	  }
	}
      }
    }
  }
  else {
    for (i=0;i<AP.NBONH;++i) {
      if (abs(AP.BH[i][0])/3==atomnumk) {
	if (strncmp(AP.IGRAPH[abs(AP.BH[i][1])/3],atomnamel,numl)==0) {
	  if (abs(AP.BH[i][1])/3 != atomnumj  && 
	      ( strncmp(AP.ITREE[abs(AP.BH[i][1])/3],"M",1)==0 || c > OMEGA || termflag != NOTERM ) ||
	      /*strncmp(AP.ITREE[abs(AP.BA[i][1])/3],"BLA",1)!=0*/ 
	      ( abs(AP.BH[i][1])/3 /*>=*/< AP.IPRES[1]-1 || abs(AP.BH[i][1])/3 /*<*/>= AP.IPRES[AP.NRES-1]-1 ) ) {
	    *atomnuml=abs(AP.BH[i][1])/3;
	    flag=S3;
	    break;
	  }
	}
      }    
      else if (abs(AP.BH[i][1])/3==atomnumk) {
	if (strncmp(AP.IGRAPH[abs(AP.BH[i][0])/3],atomnamel,numl)==0) {
	  if (abs(AP.BH[i][0])/3 != atomnumj  && 
	      ( strncmp(AP.ITREE[abs(AP.BA[i][0])/3],"M",1)==0 || c > OMEGA || termflag != NOTERM) ||
	      /*strncmp(AP.ITREE[abs(AP.BA[i][1])/3],"BLA",1)!=0*/ 
	      ( abs(AP.BH[i][1])/3 /*>=*/< AP.IPRES[1]-1 || abs(AP.BH[i][1])/3 /*<*/>= AP.IPRES[AP.NRES-1]-1 ) ) {
	    *atomnuml=abs(AP.BH[i][0])/3;
	    flag=S3;
	    break;
	  }
	}
      }
    }
  }

  if (flag!=S3) *mesg=0;

  return c;
}

int judge_dihedtype(char *atomnamej, char *atomnamek/*, char *treej, char *treek*/, char *resnamej,  char *resnamek,
		    char *atomnamei, char *atomnamel/*, char *treei, char *treel*/,int termflag){
  int i,j,k;
  int c;
  int numj=0,numk=0;
  int num1=0,num2=0;

  AADatatable *elemnt;
  AAdata dataofthisresidue;

  for (i=0;i<(sizeof(atomnamej)/sizeof(char));++i) {
    if (atomnamej[i]=='\0') break;
  }
  numj=i;

  for (i=0;i<(sizeof(atomnamek)/sizeof(char));++i) {
    if (atomnamek[i]=='\0') break;
  }
  numk=i;

  if (termflag==NTERM) {
    if ((strncmp(atomnamej,"N",1)==0 && numj==1 && strncmp(atomnamek,"CA",2)==0 && numk==2)) {
      c=PSI;
      memcpy(atomnamei,"H1",sizeof("H1"));
      memcpy(atomnamel,"C",sizeof("C"));
      return c;
    }
    else if ((strncmp(atomnamej,"CA",2)==0 && numj==2 && strncmp(atomnamek,"N",1)==0 && numk==1)) {
      c=PSI;
      memcpy(atomnamei,"C",sizeof("C"));
      memcpy(atomnamel,"H1",sizeof("H1"));
      return c;
    }
  }

  if (termflag==CTERM) {
    if ((strncmp(atomnamej,"CA",2)==0 && numj==2 && strncmp(atomnamek,"C",1)==0 && numk==1)) {
      c=PHI;
      memcpy(atomnamei,"N",sizeof("N"));
      memcpy(atomnamel,"O",sizeof("O"));
      return c;
    }
    else if ((strncmp(atomnamej,"C",1)==0 && numj==1 &&  strncmp(atomnamek,"CA",2)==0 && numk==2)) {
      c=PHI;
      memcpy(atomnamei,"O",sizeof("O"));
      memcpy(atomnamel,"N",sizeof("N"));
      return c;
    }
  }

  if ((strncmp(atomnamej,"CA",2)==0 && numj==2 && strncmp(atomnamek,"C",1)==0 && numk==1)) {
    c=PHI;
    memcpy(atomnamei,"N",sizeof("N"));
    memcpy(atomnamel,"N",sizeof("N"));
    /**********************************/
    /* memcpy(treei,"M",sizeof("M")); */
    /* memcpy(treel,"M",sizeof("M")); */
    /**********************************/
    return c;
  }
  else if ((strncmp(atomnamej,"C",1)==0 && numj==1 &&  strncmp(atomnamek,"CA",2)==0 && numk==2)) {
    c=PHI;
    memcpy(atomnamei,"N",sizeof("N"));
    memcpy(atomnamel,"N",sizeof("N"));
    /**********************************/
    /* memcpy(treei,"M",sizeof("M")); */
    /* memcpy(treel,"M",sizeof("M")); */
    /**********************************/
    return c;
  }
  else if ((strncmp(atomnamej,"N",1)==0 && numj==1 && strncmp(atomnamek,"CA",2)==0 && numk==2)) {
    c=PSI;
    memcpy(atomnamei,"C",sizeof("C"));
    memcpy(atomnamel,"C",sizeof("C"));
    /**********************************/
    /* memcpy(treei,"M",sizeof("M")); */
    /* memcpy(treel,"M",sizeof("M")); */
    /**********************************/
    return c;
  }
  else if ((strncmp(atomnamej,"CA",2)==0 && numj==2 && strncmp(atomnamek,"N",1)==0 && numk==1)) {
    c=PSI;
    memcpy(atomnamei,"C",sizeof("C"));
    memcpy(atomnamel,"C",sizeof("C"));
    /**********************************/
    /* memcpy(treei,"M",sizeof("M")); */
    /* memcpy(treel,"M",sizeof("M")); */
    /**********************************/
    return c;
  }
  else if ((strncmp(atomnamej,"N",1)==0 && numj==1 && strncmp(atomnamek,"C",1)==0 && numk==1)) {
    c=OMEGA;
    //    memcpy(atomnamei,"C",sizeof("C"));
    memcpy(atomnamei,"CA",sizeof("CA"));
    memcpy(atomnamel,"CA",sizeof("CA"));
    /**********************************/
    /* memcpy(treei,"M",sizeof("M")); */
    /* memcpy(treel,"M",sizeof("M")); */
    /**********************************/
    return c;
  }
  else if ((strncmp(atomnamej,"C",1)==0 && numj==1 && strncmp(atomnamek,"N",1)==0 && numk==1)) {
    c=OMEGA;
    memcpy(atomnamei,"CA",sizeof("CA"));
    //    memcpy(atomnamel,"C",sizeof("C"));
    memcpy(atomnamel,"CA",sizeof("CA"));
    /**********************************/
    /* memcpy(treei,"M",sizeof("M")); */
    /* memcpy(treel,"M",sizeof("M")); */
    /**********************************/
    return c;
  }

  if (strncmp(resnamej,resnamek,4)!=0) {
    printf("error: %s is not registered\n",resnamej);
    exit(1);
  }

  if((elemnt = LURAAData(resnamej,1,dataofthisresidue))==NULL) {
    printf("error: %s is not registered\n",resnamej);
    exit(1);
  }


  for (i=0;i<elemnt[0].AAdataontable.numkai;++i) {
    for (j=0;j<(sizeof(elemnt[0].AAdataontable.atomnamepair_kai[i][1])/sizeof(char));++j) 
      if (elemnt[0].AAdataontable.atomnamepair_kai[i][1][j]=='\0') break;  
    num1=j;
    for (j=0;i<(sizeof(elemnt[0].AAdataontable.atomnamepair_kai[i][2])/sizeof(char));++j) 
      if (elemnt[0].AAdataontable.atomnamepair_kai[i][2][j]=='\0') break;  
    num2=j;

    if ( strncmp(atomnamej,
		 elemnt[0].AAdataontable.atomnamepair_kai[i][1],
		 num1)==0 
	 && strncmp(atomnamek,
		    elemnt[0].AAdataontable.atomnamepair_kai[i][2],
		    num2)==0
	 ) {
      c=OMEGA+*(elemnt[0].AAdataontable.atomnamepair_kai[i][4])-'a'+1;
      memcpy(atomnamei,elemnt[0].AAdataontable.atomnamepair_kai[i][0],
	     sizeof(elemnt[0].AAdataontable.atomnamepair_kai[i][0]));
      memcpy(atomnamel,elemnt[0].AAdataontable.atomnamepair_kai[i][3],
	     sizeof(elemnt[0].AAdataontable.atomnamepair_kai[i][3]));
      return c;
    }
    else if ( strncmp(atomnamej,
		      elemnt[0].AAdataontable.atomnamepair_kai[i][2],
		      num2)==0
	      && strncmp(atomnamek,
			 elemnt[0].AAdataontable.atomnamepair_kai[i][1],
			 num1 )==0
	      ){
      c=OMEGA+*(elemnt[0].AAdataontable.atomnamepair_kai[i][4])-'a'+1;
      memcpy(atomnamei,elemnt[0].AAdataontable.atomnamepair_kai[i][3],
	     sizeof(elemnt[0].AAdataontable.atomnamepair_kai[i][3]));
      memcpy(atomnamel,elemnt[0].AAdataontable.atomnamepair_kai[i][0],
	     sizeof(elemnt[0].AAdataontable.atomnamepair_kai[i][0]));
      return c;
    }
  }

  return -1;
}
