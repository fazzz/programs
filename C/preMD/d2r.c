
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <ctype.h>

//#include "ABAb.h"
#include "ABA_hosoku.h"

#include "EF.h"
#include "PTL.h"

#include "d2r.h"

int *readd2rinput(FILE *inputfile, int *numres, int specifymove) {
  int i,j,k,d;
  int c;

  int flagsof,flagnuma,flagnum;

  int nums,n;
  int numf;

  int *resid;

  int *residdummy;
  int numresdummy;

  residdummy=(int *)gcemalloc(sizeof(int)*1);
  numresdummy=0;
  flagnum=OFF;
  flagsof=start;
  nums=0;
  numf=0;
    
  while ((c=getc(inputfile))!=-1){
    if (c=='\n') {
      if (flagnum==ON) {
	flagnum=OFF;
	if (flagsof==start) {
	  residdummy=(int *)gcerealloc(residdummy,sizeof(int)*(numresdummy+1));
	  residdummy[numresdummy]=nums;
	  nums=0;
	  ++numresdummy;
	}
	else {
	  flagsof=start;
	  j=numresdummy;
	  numresdummy+=numf-nums+1;
	  residdummy=(int *)gcerealloc(residdummy,sizeof(int)*(numresdummy+1));
	  k=0;
	  for (;j<numresdummy;++j) {
	    residdummy[j]=nums+k;
	    ++k;
	  }
	  nums=0;
	  numf=0;
	}
      }
    }
    if (c!='#') {
      if (c==' ') {
	d=0;
	n=-1;
      }
      else if (isdigit(c)) {
	flagnum=ON;
	d=(c-'0');
	++n;
	if (flagsof==start) {
	  nums=nums*10;
	  nums+=d;
	}
	if (flagsof!=start) {
	  numf=numf*10;
	  numf+=d;
	}
      }
      if (c=='-')
	flagsof=fin;
    }
  }

  if ( specifymove==OFF ) {
    (*numres)=numresdummy;
    resid=(int *)gcemalloc(sizeof(int)*(*numres));
    for (i=0;i<(*numres);++i) {
      resid[i]=residdummy[i];
    }
  }
  else {
    (*numres)=AP.NRES-numresdummy;
    resid=(int *)gcemalloc(sizeof(int)*(*numres));
    k=0;
    for (i=0;i<AP.NRES;++i) {
      for (j=0;j<numresdummy;++j) {
	if (residdummy[j]==i+1)
	  break;
      }
      if (j==numresdummy) {
	resid[k]=i+1;
	++k;
      }
    }
  }

  return resid;
}

void writed2routput(FILE *clustfileout,CLTh *clt, int numclut){
  int i,j,k,l;

   // タンパク質中の剛体数の記述
   fprintf(clustfileout, "%2d\n", numclut);

   // 剛体の原点情報の記述
   for (i=0;i<numclut;++i) {
     fprintf(clustfileout, "%5d ", clt[i].origin_atom_a);
     if ((i+1)%10==0) fprintf(clustfileout, "\n ");
   }
   fprintf(clustfileout, "\n");

   // 剛体が終点かどうか
   for (i=0; i<numclut; ++i) {
     fprintf(clustfileout, "%5d ", clt[i].terminal);
     if ((i+1)%10==0) fprintf(clustfileout, "\n ");
   }
   fprintf(clustfileout, "\n");

   // 剛体中の原子数の記述
   for (i=0; i<numclut; ++i) {
     fprintf(clustfileout, "%5d ", clt[i].num_atom_clust);
     if ((i+1)%10==0) fprintf(clustfileout, "\n ");
   }
   fprintf(clustfileout, "\n");

   // 剛体の枝の数の記述
   for (i=0;i<numclut-1; ++i) {
     fprintf(clustfileout, "%5d ",clt[i].num_branch);
     if ((i+1)%10==0) fprintf(clustfileout, "\n ");
   }
   fprintf(clustfileout, "1 ");
   fprintf(clustfileout, "\n");

   // 剛体の自由度の記述
   for (i=0; i<numclut; ++i) {
     fprintf(clustfileout, "%5d ", 3);
     if ((i+1)%10==0) fprintf(clustfileout, "\n ");
   }
   fprintf(clustfileout, "\n");

   // 剛体の終点情報の記述
   for (i=0; i<numclut; ++i) {
     for (  j=0;j<clt[i].num_branch;++j) {
       fprintf(clustfileout, "%5d ", clt[i].terminal_atom_a[j]);
     }
     if ((i+1)%10==0) fprintf(clustfileout, "\n ");
   }
   //   fprintf(clustfileout, "%5d ", clt[i].terminal_atom_a[0]);
   fprintf(clustfileout, "\n");

   // 剛体の"親"の情報の記述
   for (i=0; i<numclut; ++i)  {
     fprintf(clustfileout, "%5d ", clt[i].nNumClutOfParent);
     if ((i+1)%10==0) fprintf(clustfileout, "\n ");
   }
   fprintf(clustfileout, "\n");

   // 剛体の"子"の情報の記述
   for (i=0; i<numclut-1; ++i) {
     if (clt[i].terminal==0)  {
       fprintf(clustfileout, "%5d ",-1);
     }
     else {
       for (  j=0;j<clt[i].num_branch;++j) {
	 fprintf(clustfileout, "%5d ", clt[i].nNumClutOfChild[j]);
       }
     }
     if ((i+1)%10==0) fprintf(clustfileout, "\n ");
   }
   fprintf(clustfileout, "%5d ",-1);
   fprintf(clustfileout, "\n");

   j = 1;
   // "計算順"の情報の記述
   for (i=0; i<numclut; ++i) {
     fprintf(clustfileout, "%5d ", j);
     ++j;
     if ((i+1)%10==0) fprintf(clustfileout, "\n ");
   }
   fprintf(clustfileout, "\n");
   
   fprintf(clustfileout, "%5d ", 0);
   fprintf(clustfileout, "\n");
}

double trans_d2r(CLTh *clt, int *numclut, int **d2r, int numd2r) {
  int i,j,k,l;
  int da,db;
  int num,ns;
  int joinflag;
  int nNumClutdummy,nNumClutLast;

  da=-1;
  db=-1;
  for (i=0;i<numd2r;++i) {
    for (j=0;j<(*numclut);++j) {
      if (clt[j].origin_atom_a==d2r[i][1]) {
	k=clt[j].nNumClutOfParent-1;
	for (l=0;l<clt[k].num_branch;++l) {
	  if (clt[k].terminal_atom_a[l]==d2r[i][0]) {
	    da=k;
	    db=j;
	    break;
	  }
	}
      }
    }
    if ( da!=-1 && db!=-1 ) {
      for (j=db;j<(*numclut);++j) {
	if ( clt[j].nNumClutOfParent == db+1)
	  clt[j].nNumClutOfParent=da+1;
	else if ( clt[j].nNumClutOfParent >  db+1)
	  clt[j].nNumClutOfParent-=1;
	for (k=0;k<clt[j].num_branch;++k) {           
	  if ( clt[j].nNumClutOfChild[k] == db+1)
	    clt[j].nNumClutOfChild[k]=da+1;
	  else if ( clt[j].nNumClutOfChild[k] > db+1)
	    clt[j].nNumClutOfChild[k]-=1;
	}
      }

      if (clt[db].join > 0) {
	num=0;
	for (nNumClutdummy=db;clt[nNumClutdummy].join!=clt[db].join-1;++nNumClutdummy) {
	  num+=clt[nNumClutdummy].num_atom_clust;
	  nNumClutLast = nNumClutdummy;
	}
	for (j=0;j<clt[da].num_branch;++j) {
	  if ( clt[da].terminal_atom_a[j]==d2r[i][0]-1 ) {          //
	    clt[da].terminal_atom_a[j]=clt[db].terminal_atom_a[0];  //
	    break;
	  }
	}
	ns=j;
	for (j=0;j<clt[da].num_branch;++j) {
	  if ( clt[da].nNumClutOfChild[j]==db+1 ) {
	    clt[da].nNumClutOfChild[j]=clt[da].nNumClutOfChild[0];
	    /**********************************************************************************************************************************/
            /* clt[da].nNumClutOfChild=(int *)gcerealloc(clt[da].nNumClutOfChild,sizeof(int)*clt[da].num_branch+clt[db].num_branch-1);	      */
            /**********************************************************************************************************************************/
	    for (k=1;k<clt[db].num_branch;++k) {
	      clt[da].nNumClutOfChild[clt[da].num_branch+k-1]=clt[db].nNumClutOfChild[k];
	    }
	    break;
	  }
	}
	if ( clt[da].nNumClutOfChild[0] == -1 )
	  clt[da].terminal=clt[db].terminal;
	clt[da].num_atom_clust+=clt[db].num_atom_clust;
	/**********************************************************************************************************************************/
        /* clt[da].terminal_atom_a=(int *)gcerealloc(clt[da].terminal_atom_a,sizeof(int)*clt[da].num_branch+clt[db].num_branch-1);	  */
        /**********************************************************************************************************************************/
	for (j=1;j<clt[db].num_branch;++j)
	  clt[da].terminal_atom_a[ns+j]=clt[db].terminal_atom_a[j];            //
	clt[da].num_branch+=clt[db].num_branch-1;
	if ( clt[db].terminal == 0 && db != (*numclut)-1 ) clt[da].num_branch-=1;

	for (j=db;j<nNumClutLast;++j) {
	  cp_clustdata(&(clt[j+1]),&(clt[j]));
	}
	--(*numclut);
	joinflag=0;
	//	for (j=0;j<(*numclut);++j) ABAb_setJoin(clt,j,joinflag);
	for (j=0;j<(*numclut);++j) ABAh_setJoin(clt,j/*,joinflag*/);
	/****************************************/
        /* clt[0].join=0;		        */
	/* for(j=1; j<(*numclut); ++j) {        */
	/*   ABA_setJoin(clt,j);	        */
	/* }				        */
        /****************************************/

      }
      // In Main Chain
      else {
	/*********************************************************************************************/
        /* for (j=0;j<clt[da].num_branch;++j) {							     */
	/*   if ( clt[da].terminal_atom_a[/\*j*\/num_branch-1]==d2r[i][0] ) {			     */
	/*     clt[da].terminal_atom_a[/\*j*\/num_branch-1]=clt[db].terminal_atom_a[0];   // 	     */
	/*     break;										     */
	/*   }											     */
	/* }											     */
        /*********************************************************************************************/
	/********************************************************************************************/
        /* ns=j;										    */
	/* for (j=0;j<clt[da].num_branch;++j) {							    */
	/*   if ( clt[da].nNumClutOfChild[j]==db+1 ) {						    */
	/*     clt[da].nNumClutOfChild[j]=clt[da].nNumClutOfChild[0];				    */
	/*     for (k=1;k<clt[db].num_branch;++k) {						    */
	/*       clt[da].nNumClutOfChild[clt[da].num_branch+k-1]=clt[db].nNumClutOfChild[k];	    */
	/*     }										    */
	/*     break;										    */
	/*   }											    */
	/* }											    */
        /********************************************************************************************/
	if ( clt[da].nNumClutOfChild[0] == -1 )
	  clt[da].terminal=clt[db].terminal;
	clt[da].num_atom_clust+=clt[db].num_atom_clust;
	/***************************************************************************************/
        /* for (j=1;j<clt[db].num_branch;++j)						       */
	/*   clt[da].terminal_atom_a[ns+j]=clt[db].terminal_atom_a[j];               //	       */
        /***************************************************************************************/
	for (j=0;j<clt[db].num_branch;++j)						       
	  clt[da].terminal_atom_a[clt[da].num_branch-1+j]=clt[db].terminal_atom_a[j];               //
	for (j=0;j<clt[db].num_branch;++j)						       
	  clt[da].nNumClutOfChild[clt[da].num_branch-1+j]=clt[db].nNumClutOfChild[j];               //
	//	num_branch_da=clt[da].num_branch;
	clt[da].num_branch+=clt[db].num_branch-1;
	if ( clt[db].terminal == 0 && db != (*numclut)-1 ) clt[da].num_branch-=1;
	for (j=db;j<(*numclut);++j) {
	  cp_clustdata(&(clt[j+1]),&(clt[j]));
	}
	--(*numclut);
	joinflag=0;
	for (j=0;j<(*numclut);++j) ABAh_setJoin(clt,j/*,joinflag*/);
	/****************************************/
        /* clt[0].join=0;		        */
	/* for(j=1; j<(*numclut); ++j) {        */
	/*   ABA_setJoin(clt,j);	        */
	/* }				        */
        /****************************************/
      }
    }
  }
}

int **res2atom_BB(int *resid, int numatom, int numres, int *numd2r) {
  int i,j,k,l;
  int Nid=-1,CAid=-1,Cid=-1;
  int numrescheck;

  int **d2r;

  (*numd2r)=0;
  d2r=(int **)gcemalloc(sizeof(int *)*1);
  d2r[0]=(int *)gcemalloc(sizeof(int)*2);
  for (i=0;i<numres;++i) {
    for (j=AP.IPRES[resid[i]-1]-1;j<AP.IPRES[resid[i]]-1;++j) {
      if (strncmp(AP.IGRAPH[j],"N\0\0\0",4)==0) {
	Nid=j+1;
      }
      else if (strncmp(AP.IGRAPH[j],"CA\0\0",4)==0) {
	numrescheck=PTL_resnum(j);
	if (strncmp(AP.LABERES[numrescheck],"PRO",3)!=0 && Nid != -1) {
	  d2r=(int **)gcerealloc(d2r,sizeof(int *)*((*numd2r)+1));
	  d2r[(*numd2r)]=(int *)gcemalloc(sizeof(int)*2);
	  d2r[(*numd2r)][0]=Nid;
	  d2r[(*numd2r)][1]=j+1;
	  ++(*numd2r);
	}
	CAid=j+1;
      }
      else if (strncmp(AP.IGRAPH[j],"C\0\0\0",4)==0 ) {
	if ( CAid != -1 ) {
	  d2r=(int **)gcerealloc(d2r,sizeof(int *)*((*numd2r)+1));
	  d2r[(*numd2r)]=(int *)gcemalloc(sizeof(int)*2);
	  d2r[(*numd2r)][0]=CAid;
	  d2r[(*numd2r)][1]=j+1;
	  ++(*numd2r);
	}
	Cid=j+1;
      }
    }
    if (j<numatom-1) {
      d2r=(int **)gcerealloc(d2r,sizeof(int *)*((*numd2r)+1));
      d2r[(*numd2r)]=(int *)gcemalloc(sizeof(int)*2);
      d2r[(*numd2r)][0]=Cid;
      d2r[(*numd2r)][1]=j+1;
      ++(*numd2r);
    }
  }

  return d2r;
}

int **res2atom_SC(int *resid, int numres, int *numd2r,CLTh *clt, int numclut) {
  int i,j,k,l;

  int **d2r;

  (*numd2r)=0;
  d2r=(int **)gcemalloc(sizeof(int *)*1);
  d2r[0]=(int *)gcemalloc(sizeof(int)*2);
  for (i=0;i<numclut;++i) {
    for (j=0;j<numres;++j) {
      if (clt[i].origin_atom_a >= AP.IPRES[resid[j]-1] && clt[i].origin_atom_a < AP.IPRES[resid[j]]-1 ) {
	if ( clt[i].terminal != 0 ) {
	  for (k=0;k<clt[i].num_branch;++k) {
	    d2r=(int **)gcerealloc(d2r,sizeof(int)*((*numd2r)+1));
	    d2r[(*numd2r)]=(int *)gcemalloc(sizeof(int)*2);
	    d2r[(*numd2r)][0]=clt[i].terminal_atom_a[k];
	    d2r[(*numd2r)][1]=clt[clt[i].nNumClutOfChild[k]-1].origin_atom_a;
	    ++(*numd2r);	      
	  }
	}
      }
    }
  }

  return d2r;
}
