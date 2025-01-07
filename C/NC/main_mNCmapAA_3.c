
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>

#include "NC.h"
#include "PTL.h"
#include "EF.h"
#include "TOPO.h"

#include "netcdf_mine.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l,d,na,ii,jj;
  int resi,resj,resk;
  int numstep,HMODE=EXC;
  double QCA;
  double pi;

  int exflag=OFF;

  double dij,dik,djk,ajik,aijk,ajki;

  double *crdref;
  int numatom,numres;
  int dummy;

  double criteria=criteria_NC;
  int numncaa,numncres;
  int  *index_natatt,**ncmap,**ncmapexc,**ncmapres;
  double atom1[3],atom2[3],atomi[3],atomj[3],atomk[3];

  double *cradii_ca;
  double vec[3];
  double length;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *progname;
  char *crdfilename,*parmtopname,*crdreffilename;
  char *outputfilename;

  FILE *crdfile,*parmtop,*crdreffile;
  FILE *outputfile;

  while((c=getopt(argc,argv,"hec:"))!=-1) {
    switch(c) {
    case 'c':
      criteria=atof(optarg);
      break;
    case 'e':
      exflag=ON;
      break;
    case 'h':
      USAGE(progname);
      exit(1);
    default:
      USAGE(progname);
      exit(1);
    }
  }

  pi=acos(-1.0);

  progname=*argv;

  argc-=optind;
  argv+=optind;

  if (argc < 3) {
    USAGE(progname);
    exit(1);
  }
  crdreffilename = *argv;
  parmtopname = *++argv;
  outputfilename = *++argv;

  parmtop=efopen(parmtopname,"r");
  readParmtopL(parmtop);
  fclose(parmtop);
  numatom=AP.NATOM;
  numres=AP.NRES;

  crdref=(double *)gcemalloc(sizeof(double)*numatom*3);

  crdreffile=efopen(crdreffilename,"r");
  io_scanconf_Amber_ini(crdreffile,numatom,crdref);
  fclose(crdreffile);

  ncmap=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) ncmap[i]=(int *)gcemalloc(sizeof(int)*numatom);
  ncmapres=(int **)gcemalloc(sizeof(int *)*numres);
  for (i=0;i<numres;++i) ncmapres[i]=(int *)gcemalloc(sizeof(int)*numres);

  //  index_natatt=make_native_contact_list_aa_3_nadjacent_2(&numncaa,&numncres,crdref,numatom,numres,criteria,ncmap,ncmapres,EXC);
  index_natatt=make_native_contact_list_aa_3_nadjacent(&numncaa,&numncres,crdref,numatom,numres,criteria,ncmap,ncmapres,EXC);
  cradii_ca=(double *)gcemalloc(sizeof(double)*numncaa);

  ncmapexc=(int **)gcemalloc(sizeof(int *)*numatom);
  for (i=0;i<numatom;++i) ncmapexc[i]=(int *)gcemalloc(sizeof(int)*numatom);
  for (i=0;i<numatom;++i) for (j=0;j<numatom;++j) ncmapexc[i][j]=0;
  
  if (exflag==ON) {
    dummy=0;
    for (i=0;i<numatom;++i) {
      for (j=i+1;j<numatom;++j) {
	if (ncmap[i][j]==0) {
	  for (k=0;k<numatom;++k) {
	    resi=PTL_resnum(i);
	    resj=PTL_resnum(j);
	    resk=PTL_resnum(k);
	    if (/*resi!=resk*/abs(resi-resk)>1) {
	      if (k!=i && k!=j && strncmp(AP.IGRAPH[k],"H",1)!=0) {
		for (l=0;l<3;++l) {
		  atomi[l]=crdref[i*3+l];
		  atomj[l]=crdref[j*3+l];
		  atomk[l]=crdref[k*3+l];
		}
		dij=len(atomi,atomj);
		dik=len(atomi,atomk);
		djk=len(atomj,atomk);
		ajik=ang(atomj,atomi,atomk);
		aijk=ang(atomi,atomj,atomk);
		ajki=ang(atomj,atomk,atomi);
		if (((dik < dij)&& (ajik < 35.0*pi/180.0)) ) {
		  ncmapexc[i][j]=-1;
		  ++dummy;
		  break;
		}
	      }
	    }
	    if (/*resj!=resk*/abs(resj-resk)>1) {
	      if (k!=i && k!=j && strncmp(AP.IGRAPH[k],"H",1)!=0) {
		for (l=0;l<3;++l) {
		  atomi[l]=crdref[i*3+l];
		  atomj[l]=crdref[j*3+l];
		  atomk[l]=crdref[k*3+l];
		}
		dij=len(atomi,atomj);
		dik=len(atomi,atomk);
		djk=len(atomj,atomk);
		ajik=ang(atomj,atomi,atomk);
		aijk=ang(atomi,atomj,atomk);
		ajki=ang(atomj,atomk,atomi);
		if (((djk < dij)&& (ajik < 35.0*pi/180.0)) ) {
		  ncmapexc[i][j]=-1;
		  ++dummy;
		  break;
		}
	      }
	    }
	  }
	}
      }
    }

    /************************************************************************/
    /* if (exflag==ON) {							  */
    /*   dummy=0;								  */
    /*   for (i=0;i<numatom;++i) {					  */
    /*     for (j=i+1;j<numatom;++j) {					  */
    /* 	for (k=j+1;k<numatom;++k) {					  */
    /* 	  if ( ncmap[i][j]==0  && ncmap[i][k]==0 ) {			  */
    /* 	    for (l=0;l<3;++l) {						  */
    /* 	      atomi[l]=crdref[i*3+l];					  */
    /* 	      atomj[l]=crdref[j*3+l];					  */
    /* 	      atomk[l]=crdref[k*3+l];					  */
    /* 	    }								  */
    /* 	    dij=len(atomi,atomj);					  */
    /* 	    dik=len(atomi,atomk);					  */
    /* 	    djk=len(atomj,atomk);					  */
    /* 	    ajik=ang(atomj,atomi,atomk);				  */
    /* 	    aijk=ang(atomi,atomj,atomk);				  */
    /* 	    ajki=ang(atomj,atomk,atomi);				  */
    /* 	    if (((dik < dij) && (ajik < 35.0*pi/180.0))) {		  */
    /* 	      ncmapexc[i][j]=-1;					  */
    /* 	      ++dummy;							  */
    /* 	      break;							  */
    /* 	    }								  */
    /* 	    else if (((dij < dik) && (ajik < 35.0*pi/180.0))) {		  */
    /* 	      ncmapexc[i][k]=-1;					  */
    /* 	      ++dummy;							  */
    /* 	      break;							  */
    /* 	    }								  */
    /* 	  }								  */
    /* 	}								  */
    /*     }								  */
    /*   }								  */
    /************************************************************************/
    
    for (i=0;i<numatom;++i) {
      for (j=i+1;j<numatom;++j) {
	if (ncmapexc[i][j]==-1 ) {
	  ncmap[i][j]=-1;
	}
      }
    }
    
  }


  for (i=0;i<numres;++i) {
    for (j=0;j<numres;++j) {
      ncmapres[i][j]=-100;
    }
  }
  
  na=0;
  for (i=0;i<numatom;++i) {
    for (j=i+1;j<numatom;++j) {
      if (ncmap[i][j]==0 ) {
	length = 0.0;
	for(k=0;k<3;++k){
	  vec[k] = crdref[i*3+k]-crdref[j*3+k];
	  length += vec[k]*vec[k];
	}
	length = sqrt(length);
	printf("%2d:%d-%d %lf ",na+1,i,j,length);
	/**********************************************************/
	/* for (k=0;k<2;++k) printf("%c",AP.IGRAPH[i][k]);	    */
	/* printf(" ");					    */
	/* for (k=0;k<2;++k) printf("%c",AP.IGRAPH[j][k]);	    */
	/**********************************************************/
	printf("%d -%d  ",PTL_resnum(i),PTL_resnum(j));
	printf("\n ");
	++na;
	resi=PTL_resnum(i)-1;
	resj=PTL_resnum(j)-1;
	if (ncmapres[resi][resj]==-100) ncmapres[resi][resj]=1;
	  else ++ncmapres[resi][resj];
      }
    }
  }

  /**************************************/
  /* printf("   ",i);		        */
  /* for (i=0;i<numres;++i) {	        */
  /*   printf("%2d ",i+1);	        */
  /* }				        */
  /* printf("\n ");		        */
  /* for (i=0;i<numres;++i) {	        */
  /*   printf("%2d ",i+1);	        */
  /*   for (j=0;j<numres;++j) {	        */
  /*     printf("%2d ",ncmapres[i][j]); */
  /*   }			        */
  /*   printf("\n ");		        */
  /* }				        */
  /**************************************/

  printf("%d\n ",na);

  outputfile=efopen(outputfilename,"w");
  for (i=0;i<numres;++i) {
    for (j=0;j<numres;++j) {
      fprintf(outputfile,"%2d ",ncmapres[i][j]);
    }
  }
  fclose(outputfile);

  return 0;
}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-h] help \n");
  printf("%s [-h] reffilename parmfilename \n",progname);
}

