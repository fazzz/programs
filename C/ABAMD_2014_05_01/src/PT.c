#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "PT.h"
#include "EF.h"


int readParmtop(FILE *parmfile/*, struct AmberParm AP*/){
  int i,j;
  char *line,dummy;
  size_t len=0;
  getline(&line,&len,parmfile);
  getline(&line,&len,parmfile);
  getline(&line,&len,parmfile);
  getline(&line,&len,parmfile);
  //  fscanf(parmfile,"%4s",&AP.ITITLE);

  for (i=0;i<2;++i){
    getline(&line,&len,parmfile);
  }
  fscanf(parmfile,"%8d",&AP.NATOM);
  fscanf(parmfile,"%8d",&AP.NTYPES);
  fscanf(parmfile,"%8d",&AP.NBONH);
  fscanf(parmfile,"%8d",&AP.MBONA);
  fscanf(parmfile,"%8d",&AP.NTHETH);
  fscanf(parmfile,"%8d",&AP.MTHETA);
  fscanf(parmfile,"%8d",&AP.NPHIH);
  fscanf(parmfile,"%8d",&AP.MPHIA);
  fscanf(parmfile,"%8d",&AP.NHPARM);
  fscanf(parmfile,"%8d",&AP.NPARM);
  fscanf(parmfile,"%8d",&AP.NNB);
  AP.NEXT = AP.NNB;
  fscanf(parmfile,"%8d",&AP.NRES);
  fscanf(parmfile,"%8d",&AP.NBONA);
  fscanf(parmfile,"%8d",&AP.NTHETA);
  fscanf(parmfile,"%8d",&AP.NPHIA);
  fscanf(parmfile,"%8d",&AP.NUMBND);
  fscanf(parmfile,"%8d",&AP.NUMANG);
  fscanf(parmfile,"%8d",&AP.NPTRA);
  fscanf(parmfile,"%8d",&AP.NATYP);
  fscanf(parmfile,"%8d",&AP.NPHB);
  fscanf(parmfile,"%8d",&AP.IFPERT);
  fscanf(parmfile,"%8d",&AP.NBPER);
  fscanf(parmfile,"%8d",&AP.NGPER);
  fscanf(parmfile,"%8d",&AP.NDPER);
  fscanf(parmfile,"%8d",&AP.MBPER);
  fscanf(parmfile,"%8d",&AP.MGPER);
  fscanf(parmfile,"%8d",&AP.MDPER);
  fscanf(parmfile,"%8d",&AP.IFBOX);
  fscanf(parmfile,"%8d",&AP.NMXPS);
  fscanf(parmfile,"%8d",&AP.IFCAP);
  fscanf(parmfile,"%8d",&AP.NEXTRA);

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    for (j=0;j<4;++j) {
      fscanf(parmfile,"%c",&dummy);
      if (dummy!='\n')
	AP.IGRAPH[i][j]=dummy;
      else
	--j;
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%lf",&AP.CHRG[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%lf",&AP.AMASS[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%8d",&AP.IAC[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%8d",&AP.NUMEX[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<(AP.NTYPES)*(AP.NTYPES);++i){
    fscanf(parmfile,"%8d",&AP.ICO[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NRES;++i){
    //    fscanf(parmfile,"%4s",&AP.LABERES[i]); // 2014-07-04
    for (j=0;j<4;++j){                           // 2014-07-04
      fscanf(parmfile,"%c",&(AP.LABERES[i][j])); // 2014-07-04
    }                                            // 2014-07-04
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NRES;++i){
    fscanf(parmfile,"%8d",&AP.IPRES[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NUMBND;++i){
    fscanf(parmfile,"%lf",&AP.RK[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NUMBND;++i){
    fscanf(parmfile,"%lf",&AP.REQ[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NUMANG;++i){
    fscanf(parmfile,"%lf",&AP.TK[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NUMANG;++i){
    fscanf(parmfile,"%lf",&AP.TEQ[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NPTRA;++i){
    fscanf(parmfile,"%lf",&AP.PK[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NPTRA;++i){
    fscanf(parmfile,"%lf",&AP.PN[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NPTRA;++i){
    fscanf(parmfile,"%lf",&AP.PHASE[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATYP;++i){
    fscanf(parmfile,"%lf",&AP.SOLTY[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<(AP.NTYPES)*(AP.NTYPES+1)/2;++i){
    fscanf(parmfile,"%lf",&AP.CN1[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<(AP.NTYPES)*(AP.NTYPES+1)/2;++i){
    fscanf(parmfile,"%lf",&AP.CN2[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NBONH;++i){
    for (j=0;j<3;++j){
      fscanf(parmfile,"%8d",&AP.BH[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NBONA;++i){
    for (j=0;j<3;++j){
      fscanf(parmfile,"%d",&AP.BA[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NTHETH;++i){
    for (j=0;j<4;++j){
      fscanf(parmfile,"%8d",&AP.TH[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NTHETA;++i){
    for (j=0;j<4;++j){
      fscanf(parmfile,"%d",&AP.TA[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NPHIH;++i){
    for (j=0;j<5;++j){
      fscanf(parmfile,"%8d",&AP.PH[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NPHIA;++i){
    for (j=0;j<5;++j){
      fscanf(parmfile,"%8d",&AP.PA[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NEXT;++i){
    fscanf(parmfile,"%8d",&AP.NATEX[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NPHB;++i){
    fscanf(parmfile,"%lf",&AP.ASOL[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NPHB;++i){
    fscanf(parmfile,"%lf",&AP.BSOL[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NPHB;++i){
    fscanf(parmfile,"%lf",&AP.HBCUT[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    // fscanf(parmfile,"%4s",&AP.ISYMBL[i]); // 2014-07-04
    for (j=0;j<4;++j) {
      fscanf(parmfile,"%c",&AP.ISYMBL[i][j]); // 2014-07-04
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    //    fscanf(parmfile,"%4s",&AP.ITREE[i]);
    for (j=0;j<4;++j) {
      fscanf(parmfile,"%c",&AP.ITREE[i][j]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%8d",&AP.JOIN[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  for (i=0;i<AP.NATOM;++i){
    fscanf(parmfile,"%8d",&AP.IROTAT[i]);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  if (AP.IFBOX > 0) {
    fscanf(parmfile,"%d",&AP.IPTRES);

    fscanf(parmfile,"%d",&AP.NSPM);

    fscanf(parmfile,"%d",&AP.NSPSOL);

    for (i=0;i<AP.NSPM;++i){
      fscanf(parmfile,"%d",&AP.NSP[i]);
    }

    fscanf(parmfile,"%lf",&AP.BETA);

    for (i=0;i<3;++i){
      //      fscanf(parmfile,"%12d",&AP.BOX[i]); // 2014-07-04
      fscanf(parmfile,"%lf",&AP.BOX[i]);
    }
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  if (AP.IFCAP > 0) {
    fscanf(parmfile,"%d",&AP.NATCAP);

    fscanf(parmfile,"%lf",&AP.CUTCAP);

    fscanf(parmfile,"%lf",&AP.XCAP);

    fscanf(parmfile,"%lf",&AP.YCAP);

    fscanf(parmfile,"%lf",&AP.ZCAP);
  }

  for (i=0;i<3;++i){
    getline(&line,&len,parmfile);
  }
  if (AP.IFPERT > 0) {
    for (i=0;i<AP.NBPER;++i){
      for (j=0;j<2;++j){
	fscanf(parmfile,"%d",&AP.BPER[i][j]);
      }
    }

    for (i=0;i<AP.NBPER*2;++i){
      fscanf(parmfile,"%d",&AP.ICBPER[i]);
    }

    for (i=0;i<AP.NGPER;++i){
      for (j=0;j<3;++j){
	fscanf(parmfile,"%d",&AP.TPER[i][j]);
      }
    }

    for (i=0;i<AP.NGPER*2;++i){
      fscanf(parmfile,"%d",&AP.ICTPER[i]);
    }

    for (i=0;i<AP.NDPER;++i){
      for (j=0;j<4;++j){
	fscanf(parmfile,"%d",&AP.PPER[i][j]);
      }
    }

    for (i=0;i<AP.NDPER*2;++i){
      fscanf(parmfile,"%d",&AP.ICPPER[i]);
    }

    for (i=0;i<AP.NRES;++i){
      //      fscanf(parmfile,"%4s",&AP.LABERES[i]); // 2014-07-04
      for (j=0;j<4;++j) {                            // 2014-07-04
	fscanf(parmfile,"%c",&AP.LABERES[i][j]);     // 2014-07-04
      }                                              // 2014-07-04
    }

    for (i=0;i<AP.NATOM;++i){
      //      fscanf(parmfile,"%4s",&AP.IGRPER[i]);  // 2014-07-04
      for (j=0;j<4;++j) {                            // 2014-07-04
	fscanf(parmfile,"%c",&AP.IGRPER[i][j]);      // 2014-07-04
      }                                              // 2014-07-04
    }

    for (i=0;i<AP.NATOM;++i){
      //      fscanf(parmfile,"%4s",&AP.ISMPER[i]);  // 2014-07-04
      for (j=0;j<4;++j) {                            // 2014-07-04
	fscanf(parmfile,"%4s",&AP.ISMPER[i][j]);     // 2014-07-04
      }                                              // 2014-07-04
    }

    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.ALMPER[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.IAPER[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%d",&AP.IACPER[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.CGPER[i]);
    }
  } 

  /* if (AP.IPOL == 1) {
    fscanf(parmfile,"%18.8e",&AP.ATPOL[i]);

    fscanf(parmfile,"%18.8e",&AP.ATPOL1[i]);
    }*/

  if (AP.NPARM == 1) {
    fscanf(parmfile,"%d",&AP.NLES_NTYP);

    for (i=0;i<AP.NATOM;++i){
      //      fscanf(parmfile,"%lf",&AP.LES_TYPE[i]); // 2014-07-04
      fscanf(parmfile,"%d",&AP.LES_TYPE[i]);          // 2014-07-04
    }

    for (i=0;i<(AP.NLES_NTYP)*(AP.NLES_NTYP);++i){
      fscanf(parmfile,"%lf",&AP.LES_FAC[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.LES_CNUM[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fscanf(parmfile,"%lf",&AP.LES_ID[i]);
    }

  }

  return 0; // 2014-07-04
}



int readdihedpairs(int **atomdihedpairs, int *num) {
  int i,j;
  int phi[4],psi[4],omega[4],ipsi[4],iomega[4],fomega[4],fphi[4],sumnum;
  int kai[6][4];
  int PHIFLAG=-1,PSIFLAG=-1,OMEGAFLAG=-1,ACEFLAG,NMEFLAG;
  int numphsi,numomega,numkai,numres;
  char RESNAME[4];

  numphsi=0;numomega=0;numkai=0;
  
  for (i=0;i<AP.NATOM;++i) {
    if (strncmp(AP.ITREE[i],"M",1)==0) {
      if (strncmp(AP.IGRAPH[i],"CH3",3)==0){
	;
      }
      else if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
	if (PHIFLAG==1) {
	  PHIFLAG=2;
	  phi[PHIFLAG]=i;
	}
	if (PSIFLAG==0) {
	  PSIFLAG=1;
	  psi[PSIFLAG]=i;
	}
	if (OMEGAFLAG==2) {
	  ++numomega;
	  atomdihedpairs[1]=(int *)gcerealloc(atomdihedpairs[1],sizeof(int)*numomega*4);
	  atomdihedpairs[1][numomega*4-4]=omega[0];
	  atomdihedpairs[1][numomega*4-3]=omega[1];
	  atomdihedpairs[1][numomega*4-2]=omega[2];
	  atomdihedpairs[1][numomega*4-1]=i;
	}
	OMEGAFLAG=0;
	omega[0]=i;
      }
      else if (strncmp(AP.IGRAPH[i],"C",1)==0){
	if (PSIFLAG==1) {
	  PSIFLAG=2;
	  psi[PSIFLAG]=i;
	}
	if (OMEGAFLAG==0) {
	OMEGAFLAG=1;
	omega[OMEGAFLAG]=i;
	}
	if (PHIFLAG==2) {
	  phi[3]=i;
	  ++numphsi;
	  atomdihedpairs[0]=(int *)gcerealloc(atomdihedpairs[0],sizeof(int)*numphsi*4);
	  atomdihedpairs[0][numphsi*4-4]=phi[0];
	  atomdihedpairs[0][numphsi*4-3]=phi[1];
	  atomdihedpairs[0][numphsi*4-2]=phi[2];
	  atomdihedpairs[0][numphsi*4-1]=i;
	}
	PHIFLAG=0;
	phi[0]=i;
      }
      else if (strncmp(AP.IGRAPH[i],"N",1)==0){
	if (PHIFLAG==0) {
	  PHIFLAG=1;
	  phi[PHIFLAG]=i;
	}
	if (OMEGAFLAG==1) {
	  OMEGAFLAG=2;
	  omega[OMEGAFLAG]=i;
	}
	if (PSIFLAG==2) {
	  psi[3]=i;
 	  ++numphsi;
	  atomdihedpairs[0]=(int *)gcerealloc(atomdihedpairs[0],sizeof(int)*numphsi*4);
	  atomdihedpairs[0][numphsi*4-4]=psi[0];
	  atomdihedpairs[0][numphsi*4-3]=psi[1];
	  atomdihedpairs[0][numphsi*4-2]=psi[2];
	  atomdihedpairs[0][numphsi*4-1]=i;
	}
	PSIFLAG=0;
	psi[0]=i;
      }
    }
    numres=joinatomtores(i,RESNAME);
    if (strncmp(RESNAME,"ALA",3)==0) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	kai[0][0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CA",2)==0){
    	kai[0][1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CB",2)==0){
    	kai[0][2]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HB1",3)==0){
	++numkai;
	atomdihedpairs[2]=(int *)gcerealloc(atomdihedpairs[2],sizeof(int)*numkai*4);
	atomdihedpairs[2][numkai*4-4]=kai[0][0];
	atomdihedpairs[2][numkai*4-3]=kai[0][1];
	atomdihedpairs[2][numkai*4-2]=kai[0][2];
	atomdihedpairs[2][numkai*4-1]=i;
      }
    }
    else if (strncmp(RESNAME,"ACE",3)==0) {
      ACEFLAG=1;
      if (strncmp(AP.IGRAPH[i],"HH31",4)==0){
    	ipsi[0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CH3",3)==0){
    	ipsi[1]=i;
	iomega[0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"C",1)==0){
    	ipsi[2]=i;
	iomega[1]=i;
      }
    }
    else if (strncmp(RESNAME,"NME",3)==0) {
      NMEFLAG=1;
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
    	fomega[0]=i;
	fphi[0]=i;
      }
      if (strncmp(AP.IGRAPH[i],"CH3",3)==0){
	atomdihedpairs[4][0]=phi[2];
	atomdihedpairs[4][1]=phi[3];
	atomdihedpairs[4][2]=fomega[0];
	atomdihedpairs[4][3]=i;
	fphi[1]=i;
      }
      if (strncmp(AP.IGRAPH[i],"HH31",4)==0){
	atomdihedpairs[4][4]=phi[3];
	atomdihedpairs[4][5]=fphi[0];
	atomdihedpairs[4][6]=fphi[1];
	atomdihedpairs[4][7]=i;
      }
    }
    /*   else if (strncmp(AP.IGRAPH[i],"N",1)==0){	    */
    /* 	;						    */
    /*   }						    */
    /*   else if (strncmp(AP.IGRAPH[i],"CB",1)==0){	    */
    /* 	;						    */
    /*   }						    */
    /*   else if (strncmp(AP.IGRAPH[i],"HA1",3)==0){	    */
    /* 	;						    */
    /*   }						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"ASP",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"GLU",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"LEU",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"ILE",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"ASN",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"GLN",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"VAL",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"SER",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"THR",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"CYX",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"CYS",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"HID",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"HIE",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"HIP",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"MET",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"PRO",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"ARG",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"LYS",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"PHE",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"TYR",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"TRP",3)==0) {	    */
    /*   ;						    */
    /* }						    */
    /********************************************************/

    if (numres==1 && ACEFLAG==1) {
      if (strncmp(AP.IGRAPH[i],"N",1)==0){
	iomega[2]=i;
	atomdihedpairs[3][0]=ipsi[0];
	atomdihedpairs[3][1]=ipsi[1];
	atomdihedpairs[3][2]=ipsi[2];
	atomdihedpairs[3][3]=i;
      }
      else if (strncmp(AP.IGRAPH[i],"C",1)==0){
	iomega[3]=i;
	atomdihedpairs[3][4]=iomega[0];
	atomdihedpairs[3][5]=iomega[1];
	atomdihedpairs[3][6]=iomega[2];
	atomdihedpairs[3][7]=i;
      }
    }
  }

  num[0]=numphsi;num[1]=numomega;num[2]=numkai;
  if (ACEFLAG==1)
    num[3]=2;
  if (NMEFLAG==1)
    num[4]=2;

  sumnum=0;
  for (i=0;i<4;++i)
    sumnum+=num[i];
  
  return sumnum;
}

int joinatomtores(int numatom, char LABERES[4]) {
  int i,j;

  for (i=0;i<AP.NRES;++i) {
    if (numatom >= AP.IPRES[i]-1 ) {
      if (i==(AP.NRES-1)){
	for (j=0;j<3;++j) {
	  LABERES[j]=AP.LABERES[i][j];
	}
	LABERES[3]='\0';
      }
      else {
	if (numatom < AP.IPRES[i+1]-1) {
	  for (j=0;j<3;++j) {
	    LABERES[j]=AP.LABERES[i][j];
	  }
	  LABERES[3]='\0';
	  break;
	}
      }
    }
  }
  
  return i;
}

int writeParmtop(FILE *parmfile){
  int i,j;
  char *line;
  size_t len=0;
  fprintf(parmfile,"%%VERSION  VERSION_STAMP = V0001.000  DATE = 12/04/10  16:15:47                  \n");
  fprintf(parmfile,"%%FLAG TITLE                                                                     \n");
  fprintf(parmfile,"%%FORMAT(20a4)                                                                   \n");

  fprintf(parmfile,"                                                                                 \n");
  fprintf(parmfile,"%%FLAG POINTERS                                                                  \n");
  fprintf(parmfile,"%%FORMAT(10I8)                                                                   \n");

  fprintf(parmfile,"%8d",AP.NATOM);
  fprintf(parmfile,"%8d",AP.NTYPES);
  fprintf(parmfile,"%8d",AP.NBONH);
  fprintf(parmfile,"%8d",AP.MBONA);
  fprintf(parmfile,"%8d",AP.NTHETH);
  fprintf(parmfile,"%8d",AP.MTHETA);
  fprintf(parmfile,"%8d",AP.NPHIH);
  fprintf(parmfile,"%8d",AP.MPHIA);
  fprintf(parmfile,"%8d",AP.NHPARM);
  fprintf(parmfile,"%8d\n",AP.NPARM);
  fprintf(parmfile,"%8d",AP.NNB);
  fprintf(parmfile,"%8d",AP.NRES);
  fprintf(parmfile,"%8d",AP.NBONA);
  fprintf(parmfile,"%8d",AP.NTHETA);
  fprintf(parmfile,"%8d",AP.NPHIA);
  fprintf(parmfile,"%8d",AP.NUMBND);
  fprintf(parmfile,"%8d",AP.NUMANG);
  fprintf(parmfile,"%8d",AP.NPTRA);
  fprintf(parmfile,"%8d",AP.NATYP);
  fprintf(parmfile,"%8d\n",AP.NPHB);
  fprintf(parmfile,"%8d",AP.IFPERT);
  fprintf(parmfile,"%8d",AP.NBPER);
  fprintf(parmfile,"%8d",AP.NGPER);
  fprintf(parmfile,"%8d",AP.NDPER);
  fprintf(parmfile,"%8d",AP.MBPER);
  fprintf(parmfile,"%8d",AP.MGPER);
  fprintf(parmfile,"%8d",AP.MDPER);
  fprintf(parmfile,"%8d",AP.IFBOX);
  fprintf(parmfile,"%8d",AP.NMXPS);
  fprintf(parmfile,"%8d\n",AP.IFCAP);
  fprintf(parmfile,"%8d\n",AP.NEXTRA);

  fprintf(parmfile,"%%FLAG ATOM_NAME                                                                 \n");
  fprintf(parmfile,"%%FORMAT(20a4)                                                                   \n");
  for (i=0;i<AP.NATOM;++i){
    for (j=0;j<4;++j)
      fprintf(parmfile,"%c",AP.IGRAPH[i][j]);
    if ((i+1)%20==0) fprintf(parmfile,"\n");
  }
  if ((i)%20!=0)
    fprintf(parmfile,"\n");
  
  fprintf(parmfile,"%%FLAG CHARGE                                                                    \n");
  fprintf(parmfile,"%%FORMAT(5E16.8)                                                                 \n");
  for (i=0;i<AP.NATOM;++i){
    fprintf(parmfile,"%16.8E",AP.CHRG[i]);
    if ((i+1)%5==0) fprintf(parmfile,"\n");
  }
 if ((i)%5!=0) fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG MASS                                                                      \n");
  fprintf(parmfile,"%%FORMAT(5E16.8)                                                                 \n");
  for (i=0;i<AP.NATOM;++i){
    fprintf(parmfile,"%16.8E",AP.AMASS[i]);
    if ((i+1)%5==0) fprintf(parmfile,"\n");
  }
  if ((i)%5!=0) fprintf(parmfile,"\n");
 
  fprintf(parmfile,"%%FLAG ATOM_TYPE_INDEX                                                           \n");
  fprintf(parmfile,"%%FORMAT(10I8)                                                                   \n");
  for (i=0;i<AP.NATOM;++i){
    fprintf(parmfile,"%8d",AP.IAC[i]);
    if ((i+1)%10==0) fprintf(parmfile,"\n");
  }
  if ((i)%10!=0) fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG NUMBER_EXCLUDED_ATOMS                                                     \n");
  fprintf(parmfile,"%%FORMAT(10I8)                                                                   \n");
  for (i=0;i<AP.NATOM;++i){
    fprintf(parmfile,"%8d",AP.NUMEX[i]);
    if ((i+1)%10==0) fprintf(parmfile,"\n");
  }
  if ((i)%10!=0)   fprintf(parmfile,"\n");
 
  fprintf(parmfile,"%%FLAG NONBONDED_PARM_INDEX                                                      \n");
  fprintf(parmfile,"%%FORMAT(10I8)                               \n");
  for (i=0;i<(AP.NTYPES)*(AP.NTYPES);++i){
    fprintf(parmfile,"%8d",AP.ICO[i]);
    if ((i+1)%10==0) fprintf(parmfile,"\n");
  }
  if ((i)%10!=0)   fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG RESIDUE_LABEL                                                             \n");
  fprintf(parmfile,"%%FORMAT(20a4)                                                                   \n");
  for (i=0;i<AP.NRES;++i){
    fprintf(parmfile,"%-4s",AP.LABERES[i]);
    if ((i+1)%10==0) fprintf(parmfile,"\n");
  }
  if ((i)%10!=0) fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG RESIDUE_POINTER                                                           \n");
  fprintf(parmfile,"%%FORMAT(10I8)                                                                   \n");
  for (i=0;i<AP.NRES;++i){
    fprintf(parmfile,"%8d",AP.IPRES[i]);
    if ((i+1)%10==0) fprintf(parmfile,"\n");
  }
  if ((i)%10!=0)   fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG BOND_FORCE_CONSTANT                                                       \n");
  fprintf(parmfile,"%%FORMAT(5E16.8)                                                                 \n");
  for (i=0;i<AP.NUMBND;++i){
    fprintf(parmfile,"%16.8E",AP.RK[i]);
    if ((i+1)%5==0) fprintf(parmfile,"\n");
  }
  if ((i)%5!=0)     fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG BOND_EQUIL_VALUE                                                          \n");
  fprintf(parmfile,"%%FORMAT(5E16.8)                                                                 \n");
  for (i=0;i<AP.NUMBND;++i){
    fprintf(parmfile,"%16.8E",AP.REQ[i]);
    if ((i+1)%5==0) fprintf(parmfile,"\n");
  } 
  if ((i)%5!=0)    fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG ANGLE_FORCE_CONSTANT                                                      \n");
  fprintf(parmfile,"%%FORMAT(5E16.8)                                                                 \n");
  for (i=0;i<AP.NUMANG;++i){
    fprintf(parmfile,"%16.8E",AP.TK[i]);
    if ((i+1)%5==0) fprintf(parmfile,"\n");
  }
  if ((i)%5!=0)     fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG ANGLE_EQUIL_VALUE                                                         \n");
  fprintf(parmfile,"%%FORMAT(5E16.8)                                                                 \n");
  for (i=0;i<AP.NUMANG;++i){
    fprintf(parmfile,"%16.8E",AP.TEQ[i]);
    if ((i+1)%5==0) fprintf(parmfile,"\n");
  }
  if ((i)%5!=0)     fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG DIHEDRAL_FORCE_CONSTANT                                                   \n");
  fprintf(parmfile,"%%FORMAT(5E16.8)                                                                 \n");
  for (i=0;i<AP.NPTRA;++i){
    fprintf(parmfile,"%16.8E",AP.PK[i]);
    if ((i+1)%5==0) fprintf(parmfile,"\n");
  }
  if ((i)%5!=0)     fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG DIHEDRAL_PERIODICITY                                                      \n");
  fprintf(parmfile,"%%FORMAT(5E16.8)                                                                 \n");
  for (i=0;i<AP.NPTRA;++i){
    fprintf(parmfile,"%16.8E",AP.PN[i]);
    if ((i+1)%5==0) fprintf(parmfile,"\n");
  }
  if ((i)%5!=0)     fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG DIHEDRAL_PHASE                                                            \n");
  fprintf(parmfile,"%%FORMAT(5E16.8)                                                                 \n");
  for (i=0;i<AP.NPTRA;++i){
    fprintf(parmfile,"%16.8E",AP.PHASE[i]);
    if ((i+1)%5==0) fprintf(parmfile,"\n");
  }
  if ((i)%5!=0)     fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG SOLTY                                                                     \n");
  fprintf(parmfile,"%%FORMAT(5E16.8)                                                                 \n");
  for (i=0;i<AP.NATYP;++i){
    fprintf(parmfile,"%16.8E",AP.SOLTY[i]);
    if ((i+1)%5==0) fprintf(parmfile,"\n");
  }
  if ((i)%5!=0)     fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG LENNARD_JONES_ACOEF                                                       \n");
  fprintf(parmfile,"%%FORMAT(5E16.8)                                                                 \n");
  for (i=0;i<(AP.NTYPES)*(AP.NTYPES+1)/2;++i){
    fprintf(parmfile,"%16.8E",AP.CN1[i]);
    if ((i+1)%5==0) fprintf(parmfile,"\n");
  }
  if ((i)%5!=0)     fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG LENNARD_JONES_BCOEF                                                       \n");
  fprintf(parmfile,"%%FORMAT(5E16.8)                                                                 \n");
  for (i=0;i<(AP.NTYPES)*(AP.NTYPES+1)/2;++i){
    fprintf(parmfile,"%16.8E",AP.CN2[i]);
    if ((i+1)%5==0) fprintf(parmfile,"\n");
  }
  if ((i)%5!=0)     fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG BONDS_INC_HYDROGEN                                                        \n");
  fprintf(parmfile,"%%FORMAT(10I8)                                                                   \n");
  for (i=0;i<AP.NBONH;++i){
    for (j=0;j<3;++j){
      fprintf(parmfile,"%8d",AP.BH[i][j]);
      if ((i*3+j+1)%10==0) fprintf(parmfile,"\n");
    }
  }
  if (((i-1)*3+2+1)%5!=0)     fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG BONDS_WITHOUT_HYDROGEN                                                    \n");
  fprintf(parmfile,"%%FORMAT(10I8)                                                                   \n");
  for (i=0;i<AP.NBONA;++i){
    for (j=0;j<3;++j){
      fprintf(parmfile,"%8d",AP.BA[i][j]);
      if ((i*3+j+1)%10==0) fprintf(parmfile,"\n");
    }
  }
  if (((i-1)*3+2+1)%10!=0)     fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG ANGLES_INC_HYDROGEN                                                       \n");
  fprintf(parmfile,"%%FORMAT(10I8)                                                                   \n");
  for (i=0;i<AP.NTHETH;++i){
    for (j=0;j<4;++j){
      fprintf(parmfile,"%8d",AP.TH[i][j]);
      if ((i*4+j+1)%10==0) fprintf(parmfile,"\n");
    }
  }
  if (((i-1)*4+3+1)%10!=0)fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG ANGLES_WITHOUT_HYDROGEN                                                   \n");
  fprintf(parmfile,"%%FORMAT(10I8)                                                                   \n");
  for (i=0;i<AP.NTHETA;++i){
    for (j=0;j<4;++j){
      fprintf(parmfile,"%8d",AP.TA[i][j]);
      if ((i*4+j+1)%10==0) fprintf(parmfile,"\n");
    }
  }
  if (((i-1)*4+3+1)%10!=0) fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG DIHEDRALS_INC_HYDROGEN                                                    \n");
  fprintf(parmfile,"%%FORMAT(10I8)                                                                   \n");
  for (i=0;i<AP.NPHIH;++i){
    for (j=0;j<5;++j){
      fprintf(parmfile,"%8d",AP.PH[i][j]);
      if ((i*5+j+1)%10==0) fprintf(parmfile,"\n");
    }
  }
  if (((i-1)*5+4+1)%10!=0) fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG DIHEDRALS_WITHOUT_HYDROGEN                                                \n");
  fprintf(parmfile,"%%FORMAT(10I8)                                                                   \n");
  for (i=0;i<AP.NPHIA;++i){
    for (j=0;j<5;++j){
      fprintf(parmfile,"%8d",AP.PA[i][j]);
      if ((i*5+j+1)%10==0) fprintf(parmfile,"\n");
    }
  }
  if (((i-1)*5+4+1)%10!=0)fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG EXCLUDED_ATOMS_LIST                                                       \n");
  fprintf(parmfile,"%%FORMAT(10I8)                                                                   \n");
  for (i=0;i<AP.NEXT;++i){
    fprintf(parmfile,"%8d",AP.NATEX[i]);
    if ((i+1)%10==0) fprintf(parmfile,"\n");
  }
  if ((i)%10!=0)   fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG HBOND_ACOEF                                                               \n");
  fprintf(parmfile,"%%FORMAT(5E16.8)                                                                 \n");
  for (i=0;i<AP.NPHB;++i){
    fprintf(parmfile,"%16.8E",AP.ASOL[i]);
      if ((i+1)%10==0) fprintf(parmfile,"\n");
  }
  if ((i)%10!=0 || i==0 )   fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG HBOND_BCOEF                                                               \n");
  fprintf(parmfile,"%%FORMAT(5E16.8)                                                                 \n");
  for (i=0;i<AP.NPHB;++i){
    fprintf(parmfile,"%16.8E",AP.BSOL[i]);
    if ((i+1)%5==0) fprintf(parmfile,"\n");
  }
  if ((i)%5!=0 || i==0 ) fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG HBCUT                                                                     \n");
  fprintf(parmfile,"%%FORMAT(5E16.8)                                                                 \n");
  for (i=0;i<AP.NPHB;++i){
    fprintf(parmfile,"%16.8E",AP.HBCUT[i]);
    if ((i+1)%5==0) fprintf(parmfile,"\n");
  }
  if ((i)%5!=0 || i==0 ) fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG AMBER_ATOM_TYPE                                                           \n");
  fprintf(parmfile,"%%FORMAT(20a4)                                                                   \n");
  for (i=0;i<AP.NATOM;++i){
    fprintf(parmfile,"%-4s",AP.ISYMBL[i]);
    if ((i+1)%20==0) fprintf(parmfile,"\n");
  }
  if ((i)%20!=0 || i==0 ) fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG TREE_CHAIN_CLASSIFICATION                                                 \n");
  fprintf(parmfile,"%%FORMAT(20a4)                                                                   \n");
  for (i=0;i<AP.NATOM;++i){
    fprintf(parmfile,"%-4s",AP.ITREE[i]);
    if ((i+1)%20==0 ) fprintf(parmfile,"\n");
  }
  if ((i)%20!=0) fprintf(parmfile,"\n");
  fprintf(parmfile,"%%FLAG JOIN_ARRAY                                                                \n");
  fprintf(parmfile,"%%FORMAT(10I8)                                                                   \n");
  for (i=0;i<AP.NATOM;++i){
    fprintf(parmfile,"%8d",AP.JOIN[i]);
    if ((i+1)%10==0) fprintf(parmfile,"\n");
  }
  if ((i)%10!=0 ) fprintf(parmfile,"\n");

  fprintf(parmfile,"%%FLAG IROTAT                                                                    \n");
  fprintf(parmfile,"%%FORMAT(10I8)                                                                   \n");
  for (i=0;i<AP.NATOM;++i){
    fprintf(parmfile,"%8d",AP.IROTAT[i]);
    if ((i+1)%10==0) fprintf(parmfile,"\n");
  }
  if ((i)%10!=0) fprintf(parmfile,"\n");

  if (AP.IFBOX > 0) {
    fprintf(parmfile,"%8d",AP.IPTRES);

    fprintf(parmfile,"%8d",AP.NSPM);

    fprintf(parmfile,"%8d",AP.NSPSOL);

    for (i=0;i<AP.NSPM;++i){
      fprintf(parmfile,"%8d",AP.NSP[i]);
    }

    fprintf(parmfile,"%16.8E",AP.BETA);

    for (i=0;i<3;++i){
      // fprintf(parmfile,"%12d",AP.BOX[i]); // 2014-07-04
      fprintf(parmfile,"%lf",AP.BOX[i]);     // 2014-07-04
    }
  }


  if (AP.IFCAP > 0) {
    fprintf(parmfile,"%8d",AP.NATCAP);

    fprintf(parmfile,"%16.8E",AP.CUTCAP);

    fprintf(parmfile,"%16.8E",AP.XCAP);

    fprintf(parmfile,"%16.8E",AP.YCAP);

    fprintf(parmfile,"%16.8E",AP.ZCAP);
  }

  if (AP.IFPERT > 0) {
    for (i=0;i<AP.NBPER;++i){
      for (j=0;j<2;++j){
	fprintf(parmfile,"%8d",AP.BPER[i][j]);
      }
    }

    for (i=0;i<AP.NBPER*2;++i){
      fprintf(parmfile,"%8d",AP.ICBPER[i]);
    }

    for (i=0;i<AP.NGPER;++i){
      for (j=0;j<3;++j){
	fprintf(parmfile,"%8d",AP.TPER[i][j]);
	if ((i+3+j)%5==0) fprintf(parmfile,"\n");
      }
    }

    for (i=0;i<AP.NGPER*2;++i){
      fprintf(parmfile,"%8d",AP.ICTPER[i]);
	if ((i+3+j)%5==0) fprintf(parmfile,"\n");
    }

    for (i=0;i<AP.NDPER;++i){
      for (j=0;j<4;++j){
	fprintf(parmfile,"%8d",AP.PPER[i][j]);
      }
    }

    for (i=0;i<AP.NDPER*2;++i){
      fprintf(parmfile,"%8d",AP.ICPPER[i]);
    }

    for (i=0;i<AP.NRES;++i){
      fprintf(parmfile,"%-4s",AP.LABERES[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fprintf(parmfile,"%-4s",AP.IGRPER[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fprintf(parmfile,"%-4s",AP.ISMPER[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fprintf(parmfile,"%16.8E",AP.ALMPER[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fprintf(parmfile,"%16.8E",AP.IAPER[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fprintf(parmfile,"%8d",AP.IACPER[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fprintf(parmfile,"%16.8E",AP.CGPER[i]);
    }
  } 

  /* if (AP.IPOL == 1) {
    fprintf(parmfile,"%18.8e",AP.ATPOL[i]);

    fprintf(parmfile,"%18.8e",AP.ATPOL1[i]);
    }*/

  if (AP.NPARM == 1) {
    fprintf(parmfile,"%8d",AP.NLES_NTYP);

    for (i=0;i<AP.NATOM;++i){
      //      fprintf(parmfile,"%16.8E",AP.LES_TYPE[i]); // 2014-07-09
      fprintf(parmfile,"%24d",AP.LES_TYPE[i]);         // 2014-07-09
    }

    for (i=0;i<(AP.NLES_NTYP)*(AP.NLES_NTYP);++i){
      fprintf(parmfile,"%16.8E",AP.LES_FAC[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fprintf(parmfile,"%16.8E",AP.LES_CNUM[i]);
    }

    for (i=0;i<AP.NATOM;++i){
      fprintf(parmfile,"%16.8E",AP.LES_ID[i]);
    }

  }

  return 0; // 2014-07-04
}


int sdvdWradii(double s){
  int i;

  for (i=0;i<(AP.NTYPES)*(AP.NTYPES+1)/2;++i){
    AP.CN1[i]=AP.CN1[i]*pow(s,12);
    AP.CN2[i]=AP.CN2[i]*pow(s,6);
  }

  return 0; // 2014-07-04
}

int sd(double s){
  int i;

  for (i=0;i<AP.NUMBND;++i)
    AP.RK[i]=AP.RK[i]*s;

  for (i=0;i<AP.NUMANG;++i)
    AP.TK[i]=AP.TK[i]*s;

  for (i=0;i<AP.NPTRA;++i)
    AP.PK[i]=AP.PK[i]*s;

  for (i=0;i<AP.NATOM;++i) {
    AP.CHRG[i]=AP.CHRG[i]*sqrt(s);
  }

  for (i=0;i<(AP.NTYPES)*(AP.NTYPES+1)/2;++i){
    AP.CN1[i]=AP.CN1[i]*s;
    AP.CN2[i]=AP.CN2[i]*s;
  }

  return 0; // 2014-07-04
}

int sd_angfc(double s){
  int i;

  for (i=0;i<AP.NUMANG;++i){
    AP.TK[i]=AP.TK[i]*s;
  }

  return 0; // 2014-07-04
}


int sd_woba(double s){
  int i;

  for (i=0;i<AP.NPTRA;++i)
    AP.PK[i]=AP.PK[i]*s;

  for (i=0;i<AP.NATOM;++i)
    AP.CHRG[i]=AP.CHRG[i]*sqrt(s);

  for (i=0;i<(AP.NTYPES)*(AP.NTYPES+1)/2;++i){
    AP.CN1[i]=AP.CN1[i]*s;
    AP.CN2[i]=AP.CN2[i]*s;
  }

  return 0; // 2014-07-04
}

int sd_wobad(double s){
  int i;

  for (i=0;i<AP.NATOM;++i)
    AP.CHRG[i]=AP.CHRG[i]*sqrt(s);

  for (i=0;i<(AP.NTYPES)*(AP.NTYPES+1)/2;++i){
    AP.CN1[i]=AP.CN1[i]*s;
    AP.CN2[i]=AP.CN2[i]*s;
  }

  return 0; // 2014-07-04
}

int countatomtype(char *name,int num) {
  int i,natomtype=0;

  for (i=0;i<AP.NATOM;++i)
    if (strncmp(AP.IGRAPH[i],name,num)==0)
      ++natomtype;

  return natomtype;
      
}

int PT_countatomtype(char *name,int num, int MODE) {
  int i,natomtype=0;

  for (i=0;i<AP.NATOM;++i)
    if ((MODE==PTYES && strncmp(AP.IGRAPH[i],name,num)==0) ||
	(MODE==PTNO  && strncmp(AP.IGRAPH[i],name,num)!=0))
      ++natomtype;

  return natomtype;
      
}


int PT_resnum(int numatom) {
  int i;

  if (numatom>AP.NATOM) return -1;

  for (i=0;i<AP.NRES;++i) {
    if (AP.IPRES[i]-1 > numatom)
      return i;
  }
  return AP.NRES;
}

int iniatomnumofres(int numres) {

  return AP.IPRES[numres];
}

int numbonds(int flag){
  if (flag==INCH) return AP.NBONH+AP.MBONA;
  else return AP.MBONA;
}

int numangles(int flag){
  if (flag==INCH) return AP.NTHETH+AP.MTHETA;
  else return AP.MTHETA;
}

int numdihedas(int flag){
  if (flag==INCH) return AP.NPHIH+AP.MPHIA;
  else return AP.MPHIA;
}

int atomnamencmp(int atomnum, char *name, int numc){
  int i;
  i=strncmp(AP.IGRAPH[atomnum],name,numc);
  return i;
}

