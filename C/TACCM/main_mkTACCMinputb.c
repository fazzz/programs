#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "ABAb.h"
#include "AminoAcid.h"
#include "mkTACCMinputb.h"

#include "d2rc.h"

#include "PTL.h"
#include "FFL.h"
#include "EF.h"

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l=0,n,*m;
  int ms;
  int numatom,num_clust;
  int flag;

  int optflag=ALL,termflag=NOTERM,seq_rangeflag=OFF;
  int specifymove=OFF;

  int numres,*resid;

  int atomnumi,atomnumj,atomnumk,atomnuml;
  int numresi,numresj,numresk,numresl;
  char *atomnamej,*atomnamek;
  char *resnamej,*resnamek;

  int ***pairs;

  int *inpindex,inpnumH,inpnumA,*indexclut;
  int *numclutparent,*terminal,*origin;
  CLTb *clt;

  int mesg;

  char *line;
  size_t len=0;

  int c;
  extern char *optarg;
  extern int optind,opterr,optopt;

  char *AAdatabasefilename,*TACCMfilename,*parmfilename,*clustfilename,*seq_range_filename;
  FILE *inputfile,*TACCMfile,*parmfile,*clustfile,*seq_range_file;

  char *progname;
  int opt_idx=1;

  struct option long_opt[] = {
    {"h",0,NULL,'h'},
    {"P",0,NULL,'P'},
    {"O",0,NULL,'O'},
    {"seq_range",1,NULL,'s'},
    {"K",1,NULL,'K'},
    {0,0,0,0}
  };

  while((c=getopt_long(argc,argv,"hPOs:K:",long_opt,&opt_idx))!=-1) {
    switch(c) {
    case 'P':
      optflag=PSI;
      break;
    case 'O':
      optflag=OMEGA;
      break;
    case 's':
      seq_range_filename=optarg;
      seq_rangeflag=ON;
      break;
    case 'K':
      optflag=atoi(optarg)+OMEGA;
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
  AAdatabasefilename  = *argv;
  clustfilename  = *++argv;
  parmfilename   = *++argv;
  TACCMfilename  = *++argv;

  if (seq_rangeflag==ON) {
    seq_range_file=efopen(seq_range_filename,"r");
    resid=readd2rcinput(seq_range_file,&numres,specifymove);
    fclose(seq_range_file);
  }

  readAADataBase(AAdatabasefilename);

  parmfile=efopen(parmfilename,"r");
  readParmtopL(parmfile);
  fclose(parmfile); 
  numatom=AP.NATOM;

  clustfile=efopen(clustfilename,"r");
  clt=ABAbp_clustscan(clustfile,&num_clust);
  fclose(clustfile);

  numclutparent=(int *)gcemalloc(sizeof(int)*num_clust);
  terminal=(int *)gcemalloc(sizeof(int)*num_clust);
  origin=(int *)gcemalloc(sizeof(int)*num_clust);
  for (i=0;i<num_clust;++i) {
    numclutparent[i]=clt[i].nNumClutOfParent;
    terminal[i]=clt[i].terminal_atom_a[0];
    origin[i]=clt[i].origin_atom_a;
  }

  indexclut=(int *)gcemalloc(sizeof(int)*(AP.NPHIH+AP.MPHIA));
  inpindex=ffL_make_inpindex(&inpnumH,&inpnumA,indexclut,num_clust,numclutparent,terminal,origin);

  atomnamej=(char *)gcemalloc(sizeof(char)*5);
  atomnamek=(char *)gcemalloc(sizeof(char)*5);
  resnamej=(char *)gcemalloc(sizeof(char)*5);
  resnamek=(char *)gcemalloc(sizeof(char)*5);

  pairs=(int ***)gcemalloc(sizeof(int **)*3); // PHIPSI,OMEGA,KAI
  for (i=0;i<3;++i) {
    pairs[i]=(int **)gcemalloc(sizeof(int *)*1);
  }
  //  for (i=0;i<3;++i) {
  //    pairs[i][0]=(int *)gcemalloc(sizeof(int)*5);
  //  }

  m=(int *)gcemalloc(sizeof(int)*3);

  //  n=1;
  for (i=0;i<3;++i) m[i]=1;
  for (i=0;i<num_clust;++i) {
    if ( clt[i].terminal!=0 ) {
      for (j=0;j<clt[i].num_branch;++j) {
	atomnumj=clt[i].terminal_atom_a[j]-1;
	atomnumk=clt[clt[i].nNumClutOfChild[j]-1].origin_atom_a-1;

	memcpy(atomnamej,AP.IGRAPH[atomnumj],sizeof(AP.IGRAPH[atomnumj]));
	memcpy(atomnamek,AP.IGRAPH[atomnumk],sizeof(AP.IGRAPH[atomnumk]));
	numresj=PTL_joinatomtores(atomnumj,resnamej);
	numresk=PTL_joinatomtores(atomnumk,resnamek);

	if ( (PTL_which_include(numresj,resid,numres) == 0 &&
	      PTL_which_include(numresk,resid,numres) == 0  ) || seq_rangeflag==OFF) {
	
	  termflag=NOTERM;
	  if (atomnumj<AP.IPRES[1]-1 ) termflag=NTERM;
	  else if ( atomnumk>(AP.IPRES[AP.NRES-1]-1)) termflag=CTERM;

	  c=mkTACCMinput_set_atomnumi_atomnuml(atomnumj,atomnumk,
					       atomnamej,atomnamek,
					       resnamej,resnamek,
					       &atomnumi,&atomnuml,&mesg,termflag);

	  if ( mesg==0  ) {
	    printf("error: clut=%d ",n);
	    if (c==PHI)
	      printf("PHI  \n");
	    else if (c==PSI)
	      printf("PSI  \n");
	    else if (c==OMEGA)
	      printf("OMEGA\n");
	    else 
	      printf("KAI%d",c-2);
	  }
	
	  if ( mesg!=0 && c <= optflag ) {
	    flag=OFF;
	    for (k=0;k<AP.NPHIH;++k) {
	      if (atomnumi==abs(AP.PH[k][0])/3 && atomnumj==abs(AP.PH[k][1])/3  &&
		  atomnumk==abs(AP.PH[k][2])/3 && atomnuml==abs(AP.PH[k][3])/3  ) {
		n=indexclut[k]-1;
		flag=ON;
		break;
	      }
	    }
	    if ( flag==OFF ) {
	      for (k=AP.NPHIH;k<AP.NPHIH+AP.MPHIA;++k) {
		if (atomnumi==abs(AP.PA[k-AP.NPHIH][0])/3 && atomnumj==abs(AP.PA[k-AP.NPHIH][1])/3  &&
		    atomnumk==abs(AP.PA[k-AP.NPHIH][2])/3 && atomnuml==abs(AP.PA[k-AP.NPHIH][3])/3  ) {
		  n=indexclut[k]-1;
		  flag=ON;
		  break;
		}
	      }
	    }

	    if (flag==OFF) printf("error: numclut=%4d\n",i);

	    //	  fprintf(TACCMfile,"%4d %4d %4d %4d %4d\n",atomnumi+1,atomnumj+1,atomnumk+1,atomnuml+1,n);
	    /*************************************************************/
            /* printf("%4d(%4s) %4d(%4s) %4d(%4s) %4d(%4s) %4d ",	 */
	    /* 	   atomnumi+1,AP.IGRAPH[atomnumi],			 */
	    /* 	   atomnumj+1,AP.IGRAPH[atomnumj],			 */
	    /* 	   atomnumk+1,AP.IGRAPH[atomnumk],			 */
	    /* 	   atomnuml+1,AP.IGRAPH[atomnuml],n);			 */
            /*************************************************************/
	    /************************************************/
            /* if (c==PHI)				    */
	    /*   printf("PHI  ");			    */
	    /* else if (c==PSI)				    */
	    /*   printf("PSI  ");			    */
	    /* else if (c==OMEGA)			    */
	    /*   printf("OMEGA");			    */
	    /* else					    */
	    /*   printf("KAI%d ",c-2);			    */
	    /* printf(" %s-%s\n",resnamej,resnamek);	    */
            /************************************************/
	    if (c==PHI || c==PSI) {
	      pairs[0]=(int **)gcerealloc(pairs[0],sizeof(int *)*m[0]);
	      pairs[0][m[0]-1]=(int *)gcemalloc(sizeof(int)*5);
	      pairs[0][m[0]-1][0]=atomnumi+1;
	      pairs[0][m[0]-1][1]=atomnumj+1;
	      pairs[0][m[0]-1][2]=atomnumk+1;
	      pairs[0][m[0]-1][3]=atomnuml+1;
	      pairs[0][m[0]-1][4]=n;
	      ++m[0];
	    }
	    else if (c==OMEGA) {
	      pairs[1]=(int **)gcerealloc(pairs[1],sizeof(int *)*m[1]);
	      pairs[1][m[1]-1]=(int *)gcemalloc(sizeof(int)*5);
	      pairs[1][m[1]-1][0]=atomnumi+1;
	      pairs[1][m[1]-1][1]=atomnumj+1;
	      pairs[1][m[1]-1][2]=atomnumk+1;
	      pairs[1][m[1]-1][3]=atomnuml+1;
	      pairs[1][m[1]-1][4]=n;
	      ++m[1];
	    }
	    else {
	      pairs[2]=(int **)gcerealloc(pairs[2],sizeof(int *)*m[2]);
	      pairs[2][m[2]-1]=(int *)gcemalloc(sizeof(int)*5);
	      pairs[2][m[2]-1][0]=atomnumi+1;
	      pairs[2][m[2]-1][1]=atomnumj+1;
	      pairs[2][m[2]-1][2]=atomnumk+1;
	      pairs[2][m[2]-1][3]=atomnuml+1;
	      pairs[2][m[2]-1][4]=n;
	      ++m[2];
	    }
	  }
	}
	//	++n;
      }
    }
  }

  TACCMfile=efopen(TACCMfilename,"w");
  ms=0; for (i=0;i<3;++i) ms+=m[i]-1;
  fprintf(TACCMfile,"%4d\n",ms);
  for (i=0;i<3;++i) 
    for (j=0;j<m[i]-1;++j) 
      fprintf(TACCMfile,"%4d %4d %4d %4d %4d\n",
	      pairs[i][j][0],pairs[i][j][1],pairs[i][j][2],pairs[i][j][3],pairs[i][j][4]);
  fclose(TACCMfile);

  return 0;

}

int USAGE(char *progname) {
  printf("USAGE:\n");
  printf("[-O] Omegaflag \n");
  printf("[-h] help \n");
  printf("%s [-h] inputfilename clustfilename parmfilename TACCMfilename \n",progname);
}
