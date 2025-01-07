
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "PT.h"
#include "EF.h"

int readdihedpairs(int **atomdihedpairs, int numphsi, int numomega, int numkai) {
  int i,j;
  int phi[4],psi[4],omega[4],ipsi[4],finphi[4];
  int kai[6][4];
  int PHIFLAG=0,PSIFLAG=0,OMEGAFLAG=0;
  
  for (i=0;i<AP.NATOM;++i) {
    if (strncmp(AP.ITREE[i],"M",1)==0) {
      if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
	if (PHIFLAG==1) {
	  PHIFLAG=2;
	  phi[PHIFLAG]=i;
	}
	if (PSIFLAG==0) {
	  PSIFLAG=1;
	  psi[PSIFLAG]=i;
	}
	if (OMEGAFLAG==2) {
	  atomdihedpairs[1]=(int *)gcerealloc(atomdihedpairs[1],sizeof(int)*numomega*4);
	  atomdihedpairs[1][numomega-4]=omega[0];
	  atomdihedpairs[1][numomega-3]=omega[1];
	  atomdihedpairs[1][numomega-2]=omega[2];
	  atomdihedpairs[1][numomega-1]=i;
	}
	++numomega;
	OFLAG=0;
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
	  atomdihedpairs[0]=(int *)gcerealloc(atomdihedpairs[0],sizeof(int)*numphsi*4);
	  atomdihedpairs[0][numphsi-4]=phi[0];
	  atomdihedpairs[0][numphsi-3]=phi[1];
	  atomdihedpairs[0][numphsi-2]=phi[2];
	  atomdihedpairs[0][numphsi-1]=i;
	}
	++numphsi;
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
	  atomdihedpairs[0]=(int *)gcerealloc(atomdihedpairs[0],sizeof(int)*numphsi*4);
	  atomdihedpairs[0][numphsi-4]=psi[0];
	  atomdihedpairs[0][numphsi-3]=psi[1];
	  atomdihedpairs[0][numphsi-2]=psi[2];
	  atomdihedpairs[0][numphsi-1]=i;
	}
	++numphsi;
	PSIFLAG=0;
	psi[0]=i;
      }
    }
    /********************************************************/
    /* if (strncmp(AP.RESNAME[i],"GLY",3)==0) {		    */
    /*   ;						    */
    /* }						    */
    /* else if (strncmp(AP.RESNAME[i],"ALA",3)==0) {	    */
    /*   if (strncmp(AP.IGRAPH[i],"CA",1)==0){		    */
    /* 	;						    */
    /*   }						    */
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
    
  }
  

  
  return 0;
}


