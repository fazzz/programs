#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <netcdf.h>
#include <ctype.h>
#include <getopt.h>

#include "ABA.h"
#include "preABAMD.h"

#include "PTL.h"
#include "EF.h"
#include "gc.h"

#define ON 1
#define OFF 0

#define ATOM 0
#define RESIDUE 1

#define start 0
#define fin   1

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,m,o,ns,d;
  int numclut,nNumClutdummy,num,nNumClutLast;
  int flag,joinflag,flagsof,flagnum,flagINC;

  int nums,n;
  int numf;

  int typeflag=ATOM;
  int specifymove=OFF;
  int sidechainmove=ON;

  int Nid,CAid,Cid;

  int numd2r;
  int numatom;
  int numres;
  int numresdummy;
  int *d2r;
  int *d2r_c;
  int da,db;
  int *resid;
  int *residdummy;

  CLT *clt;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *inputfilename,*clustfileinname,*parmfilename;
  char *clustfileoutname;

  FILE *inputfile,*clustfilein,*parmfile;
  FILE *clustfileout;

  char *progname;

  int opt_idx=1;

  struct option long_opt[] = {
    {"res",0,NULL,'e'},
    //    {"backbone",0,NULL,'b'},
    {"sidechain",0,NULL,'s'},
    {"move",0,NULL,'m'},
    {"h",0,NULL,'h'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hesm",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'e':
      typeflag=RESIDUE;
      break;
    case 'm':
      specifymove=ON;
      break;
    case 's':
      sidechainmove=OFF;
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 4) {
    USAGE(progname);
    exit(1);
  }
  inputfilename    = *argv;
  parmfilename  = *++argv;
  clustfileinname  = *++argv;
  clustfileoutname = *++argv;

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;

  clustfilein=efopen(clustfileinname,"r");
  clt=ABAp_clustscan(clustfilein,&numclut);
  fclose(clustfilein);
  
  joinflag=0;
  for (i=0;i<numclut;++i) ABA_setJoin(clt,i,joinflag);

  inputfile=efopen(inputfilename,"r");
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
    numres=numresdummy;
    resid=(int *)gcemalloc(sizeof(int)*numres);
    for (i=0;i<numres;++i) {
      resid[i]=residdummy[i];
    }
  }
  else {
    numres=AP.NRES-numresdummy;
    resid=(int *)gcemalloc(sizeof(int)*numres);
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

  if (sidechainmove==ON) {
    numd2r=0;
    d2r=(int *)gcemalloc(sizeof(int)*1);
    for (i=0;i<numres;++i) {
      for (j=AP.IPRES[resid[i]-1]-1;j<AP.IPRES[resid[i]]-1;++j) {
	if (strncmp(AP.IGRAPH[j],"N\0\0\0",4)==0) {
	  Nid=j+1;
	}
	else if (strncmp(AP.IGRAPH[j],"CA\0\0",4)==0) {
	  flagINC=OFF;
	  for (k=0;k<numclut;++k) {
	    if (clt[k].origin_atom_a=Nid ) {
	      if ( clt[k].terminal != 0) {
		for (l=0;l<clt[k].num_branch;++l) {
		  if ( clt[clt[k].nNumClutOfChild[l]-1].origin_atom_a==j+1 ) 
		    flagINC=ON;
		  for (m=0;m<clt[clt[k].nNumClutOfChild[l]-1].num_branch;++m) {
		    if ( clt[clt[k].nNumClutOfChild[l]-1].terminal_atom_a[m]==j+1 ) 
		      flagINC=ON;
		  }
		}
	      }
	    }
	  }
	  for (k=0;k<numclut;++k) {
	    for (o=0;n<clt[k].num_branch;++o) {
	      if (clt[k].terminal_atom_a[o]=Nid ) {
		if ( clt[k].terminal != 0) {
		  for (l=0;l<clt[k].num_branch;++l) {
		    if ( clt[clt[k].nNumClutOfChild[l]-1].origin_atom_a==j+1 ) 
		      flagINC=ON;
		    for (m=0;m<clt[clt[k].nNumClutOfChild[l]-1].num_branch;++m) {
		      if ( clt[clt[k].nNumClutOfChild[l]-1].terminal_atom_a[m]==j+1 ) 
			flagINC=ON;
		    }
		  }
		}
	      }
	    }
	  }
	  if (flagINC==ON) {
	    d2r=(int *)gcerealloc(d2r,sizeof(int)*(numd2r+1)*2);
	    d2r[numd2r*2+0]=Nid;
	    d2r[numd2r*2+1]=j+1;
	    ++numd2r;
	  }
	  CAid=j+1;
	}
	else if (strncmp(AP.IGRAPH[j],"C\0\0\0",4)==0) {
	  flagINC=OFF;
	  for (k=0;k<numclut;++k) {
	    if (clt[k].origin_atom_a=CAid ) {
	      if ( clt[k].terminal != 0) {
		for (l=0;l<clt[k].num_branch;++l) {
		  if ( clt[clt[k].nNumClutOfChild[l]-1].origin_atom_a==j+1 ) 
		    flagINC=ON;
		  for (m=0;m<clt[clt[k].nNumClutOfChild[l]-1].num_branch;++m) {
		    if ( clt[clt[k].nNumClutOfChild[l]-1].terminal_atom_a[m]==j+1 ) 
		      flagINC=ON;
		  }
		}
	      }
	    }
	  }
	  for (k=0;k<numclut;++k) {
	    for (o=0;o<clt[k].num_branch;++o) {
	      if (clt[k].terminal_atom_a[o]=CAid ) {
		if ( clt[k].terminal != 0) {
		  for (l=0;l<clt[k].num_branch;++l) {
		    if ( clt[clt[k].nNumClutOfChild[l]-1].origin_atom_a==j+1 ) 
		      flagINC=ON;
		    for (m=0;m<clt[clt[k].nNumClutOfChild[l]-1].num_branch;++m) {
		      if ( clt[clt[k].nNumClutOfChild[l]-1].terminal_atom_a[m]==j+1 ) 
			flagINC=ON;
		    }
		  }
		}
	      }
	    }
	  }
	  if (flagINC==ON) {
	    d2r=(int *)gcerealloc(d2r,sizeof(int)*(numd2r+1)*2);
	    d2r[numd2r*2+0]=CAid;
	    d2r[numd2r*2+1]=j+1;
	    ++numd2r;
	  }
	  Cid=j+1;
	}
      }
      if (j<numatom-1) {
	flagINC=OFF;
	for (k=0;k<numclut;++k) {
	  if (clt[k].origin_atom_a=Cid ) {
	    if ( clt[k].terminal != 0) {
	      for (l=0;l<clt[k].num_branch;++l) {
		if ( clt[clt[k].nNumClutOfChild[l]-1].origin_atom_a==j+1 ) 
		  flagINC=ON;
		for (m=0;m<clt[clt[k].nNumClutOfChild[l]-1].num_branch;++m) {
		  if ( clt[clt[k].nNumClutOfChild[l]-1].terminal_atom_a[m]==j+1 ) 
		    flagINC=ON;
		}
	      }
	    }
	  }
	}
	for (k=0;k<numclut;++k) {
	  for (o=0;o<clt[k].num_branch;++o) {
	    if (clt[k].terminal_atom_a[o]=Cid ) {
	      if ( clt[k].terminal != 0) {
		for (l=0;l<clt[k].num_branch;++l) {
		  if ( clt[clt[k].nNumClutOfChild[l]-1].origin_atom_a==j+1 ) 
		    flagINC=ON;
		  for (m=0;m<clt[clt[k].nNumClutOfChild[l]-1].num_branch;++m) {
		    if ( clt[clt[k].nNumClutOfChild[l]-1].terminal_atom_a[m]==j+1 ) 
		      flagINC=ON;
		  }
		}
	      }
	    }
	  }
	}
	if (flagINC==ON) {
	  d2r=(int *)gcerealloc(d2r,sizeof(int)*(numd2r+1)*2);
	  d2r[numd2r*2+0]=Cid;
	  d2r[numd2r*2+1]=j+1;
	  ++numd2r;
	}
      }
    }
  }
  else {
    numd2r=0;
    d2r=(int *)gcemalloc(sizeof(int)*2);
    for (i=0;i<numclut;++i) {
      for (j=0;j<numres;++j) {
	if (clt[i].origin_atom_a >= AP.IPRES[resid[j]-1] && clt[i].origin_atom_a < AP.IPRES[resid[j]]-1 ) {
	  if ( clt[i].terminal != 0 ) {
	    for (k=0;k<clt[i].num_branch;++k) {
	      d2r=(int *)gcerealloc(d2r,sizeof(int)*(numd2r+1)*2);
	      d2r[numd2r*2+0]=clt[i].terminal_atom_a[k];
	      d2r[numd2r*2+1]=clt[clt[i].nNumClutOfChild[k]-1].origin_atom_a;
	      ++numd2r;	      
	    }
	  }
	}
      }
    }
  }
  fclose(inputfile);

  da=-1;
  db=-1;
  for (i=0;i<numd2r;++i) {
    for (j=0;j<numclut;++j) {
      if (clt[j].origin_atom_a==d2r[i*2+1]) {
	k=clt[j].nNumClutOfParent-1;
	for (l=0;l<clt[k].num_branch;++l) {
	  if (clt[k].terminal_atom_a[l]==d2r[i*2+0]) {
	    da=k;
	    db=j;
	    break;
	  }
	}
      }
    }
    for (j=db;j<numclut;++j) {
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

    if ( da==-1 || db==-1 ) {
      printf("error\n");
      exit(1);
    }

    // In Side Chain
    if (clt[db].join > 0) {
      num=0;
      for (nNumClutdummy=db;clt[nNumClutdummy].join!=clt[db].join-1;++nNumClutdummy) {
	num+=clt[nNumClutdummy].num_atom_clust;
	nNumClutLast = nNumClutdummy;
      }
      for (j=0;j<clt[da].num_branch;++j) {
	if ( clt[da].terminal_atom_a==d2r[i*2+0]-1 ) {
	  clt[da].terminal_atom_a[j]=clt[db].terminal_atom_a[0];
	  break;
	}
      }
      ns=j;
      for (j=0;j<clt[da].num_branch;++j) {
	if ( clt[da].nNumClutOfChild[j]==db+1 ) {
	  clt[da].nNumClutOfChild[j]=clt[da].nNumClutOfChild[0];
	  for (k=1;k<clt[db].num_branch;++k) {
	    clt[da].nNumClutOfChild[clt[da].num_branch+k-1]=clt[db].nNumClutOfChild[k];
	  }
	  break;
	}
      }
      if ( clt[da].nNumClutOfChild[0] == -1 )
	clt[da].terminal=clt[db].terminal;
      clt[da].num_atom_clust+=clt[db].num_atom_clust;
      for (j=1;j<clt[db].num_branch;++j)
	clt[da].terminal_atom_a[ns+j]=clt[db].terminal_atom_a[j];
      clt[da].num_branch+=clt[db].num_branch-1;
      if ( clt[db].terminal == 0 && db != numclut-1 ) clt[da].num_branch-=1;

      for (j=db;j<nNumClutLast;++j) {
	cp_clustdata(&(clt[j+1]),&(clt[j]));
	//	if (j!=db)
	/**************************************************/
        /* clt[j].nNumClutOfParent-=1;			  */
	/* for (k=0;k<clt[j].num_branch;++k)		  */
	/*   if ( clt[j].nNumClutOfChild[k] != -1)	  */
	/*     clt[j].nNumClutOfChild[k]-=1;		  */
        /**************************************************/
      }
      /****************************************************/
      /* for (j=nNumClutLast;j<numclut;++j) {		  */
      /* 	clt[j].nNumClutOfParent-=1;		  */
      /* 	for (k=0;k<clt[j].num_branch;++k)	  */
      /* 	  if ( clt[j].nNumClutOfChild[k] != -1)	  */
      /* 	    clt[j].nNumClutOfChild[k]-=1;	  */
      /* }						  */
      /****************************************************/
      --numclut;
      joinflag=0;
      for (j=0;j<numclut;++j) ABA_setJoin(clt,j,joinflag);
    }
    // In Main Chain
    else {
      for (j=0;j<clt[da].num_branch;++j) {
	if ( clt[da].terminal_atom_a[j]==d2r[i*2+0] ) {
	  clt[da].terminal_atom_a[j]=clt[db].terminal_atom_a[0];
	  break;
	}
      }
      ns=j;
      for (j=0;j<clt[da].num_branch;++j) {
	if ( clt[da].nNumClutOfChild[j]==db+1 ) {
	  clt[da].nNumClutOfChild[j]=clt[da].nNumClutOfChild[0];
	  for (k=1;k<clt[db].num_branch;++k) {
	    clt[da].nNumClutOfChild[clt[da].num_branch+k-1]=clt[db].nNumClutOfChild[k];
	  }
	  break;
	}
      }
      if ( clt[da].nNumClutOfChild[0] == -1 )
	clt[da].terminal=clt[db].terminal;
      clt[da].num_atom_clust+=clt[db].num_atom_clust;
      for (j=1;j<clt[db].num_branch;++j)
	clt[da].terminal_atom_a[ns+j]=clt[db].terminal_atom_a[j];
      clt[da].num_branch+=clt[db].num_branch-1;
      if ( clt[db].terminal == 0 && db != numclut-1 ) clt[da].num_branch-=1;
      for (j=db;j<numclut;++j) {
	cp_clustdata(&(clt[j+1]),&(clt[j]));
	//	if (j!=db)
	/**************************************************/
        /* clt[j].nNumClutOfParent-=1;			  */
	/* for (k=0;k<clt[j].num_branch;++k)		  */
	/*   if ( clt[j].nNumClutOfChild[k] != -1)	  */
	/*     clt[j].nNumClutOfChild[k]-=1;		  */
        /**************************************************/
      }
      --numclut;
      joinflag=0;
      for (j=0;j<numclut;++j) ABA_setJoin(clt,j,joinflag);
    }
  }

  clustfileout=efopen(clustfileoutname,"w");
  fprintf(clustfileout,"%4d ",numclut);
  fprintf(clustfileout,"\n");
  for(i=0;i<numclut;++i) {
    fprintf(clustfileout,"%4d ",clt[i].origin_atom_a);
    for(j=1;j<clt[i].num_branch;++j) fprintf(clustfileout,"     ");
    if (i%10 == 0 && i>0)
      fprintf(clustfileout,"\n");
  }
  fprintf(clustfileout,"\n");
  for(i=0;i<numclut;++i) {
    fprintf(clustfileout,"%4d ",clt[i].terminal);
    for(j=1;j<clt[i].num_branch;++j) fprintf(clustfileout,"     ");
    if (i%10 == 0 && i>0)
      fprintf(clustfileout,"\n");
  }
  fprintf(clustfileout,"\n");
  for(i=0;i<numclut;++i) {
    fprintf(clustfileout,"%4d ",clt[i].num_atom_clust);
    for(j=1;j<clt[i].num_branch;++j) fprintf(clustfileout,"     ");
    if (i%10 == 0 && i>0)
      fprintf(clustfileout,"\n");
  }
  fprintf(clustfileout,"\n");
  for(i=0;i<numclut;++i) {
    fprintf(clustfileout,"%4d ",clt[i].num_branch);
    for(j=1;j<clt[i].num_branch;++j) fprintf(clustfileout,"     ");
    if (i%10 == 0 && i>0)
      fprintf(clustfileout,"\n");
  }
  fprintf(clustfileout,"\n");
  for(i=0;i<numclut;++i) {
    fprintf(clustfileout,"%4d ",3);
    for(j=1;j<clt[i].num_branch;++j) fprintf(clustfileout,"     ");
    if (i%10 == 0 && i>0)
      fprintf(clustfileout,"\n");
  }
  fprintf(clustfileout,"\n");
  for(i=0;i<numclut;++i) {
    for(j=0;j<clt[i].num_branch;++j) fprintf(clustfileout,"%4d ",clt[i].terminal_atom_a[j]);
    if (i%10 == 0 && i>0)
      fprintf(clustfileout,"\n");
  }
  fprintf(clustfileout,"\n");
  for (i=0;i<numclut;++i) {
    fprintf(clustfileout, "%4d ", clt[i].nNumClutOfParent);
    for(j=1;j<clt[i].num_branch;++j) fprintf(clustfileout,"     ");
    if (i%10 == 0 && i>0)
      fprintf(clustfileout,"\n");
  }
  fprintf(clustfileout,"\n");
  for (i=0;i<numclut;++i) {
    for(j=0;j<clt[i].num_branch;++j) fprintf(clustfileout, "%4d ", clt[i].nNumClutOfChild[j]);
    if (i%10 == 0 && i>0)
      fprintf(clustfileout,"\n");
  }
  fprintf(clustfileout,"\n");
  for (i=0;i<numclut;++i) {
    fprintf(clustfileout, "%4d ", i+1);
    for(j=1;j<clt[i].num_branch;++j) fprintf(clustfileout,"     ");
    if (i%10 == 0 && i>0)
      fprintf(clustfileout,"\n");
  }
  fprintf(clustfileout,"\n");
  fprintf(clustfileout,"0 \n");
  fclose(clustfileout);

  return 0;
  
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfileinname clustfileoutputfilename \n",progname);
}
