#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "EF.h"
#include "PT.h"
#include "RBD.h"
#include "LA.h"

int main(int argc, char *argv[]) {
  int i,j,k;

  int natom1,natom2,natom3,natomi,natomj,numatom,numbody,nparent,numatom_body;
  double *coord,*coord_body,*velo,*crd_org,*Tmat;
  double *mass,*Inertia,*InvInertia;
  double *velo_lab,*velo_body,*velo_temp;
  double *coord_body_temp;
  double *mom,*cv,*og;
  double *Rmat,*Rotmat,*Rotmati,*Rotmatj;
  double *Transmat,*Transmati_j,*TTransmati_j;
  double *spvelo,*spvelo_pre,*Tspvelo;
  double *dtheta;
  double *theta,*TRmat,*TRtheta;
  double *omeg_term,*velo_term,*velo_com;
  double det;
  RigidBodyColl RBC;

  char *inputfilename1,*inputfilename2,*inputfilename3,*outputfilename;
  FILE *inputfile1,*inputfile2,*inputfile3, *outputfile;

  char *line;
  size_t len=0;
  
  if (argc < 5) {
    printf("USAGE: tra_rst_cMDtotaMD inputfilename1(rstcMD) inputfilename2(clt) inputfilename3(parm) outputfilename(rsttaMD)\n");
    exit(1);
  }
  inputfilename1 = *++argv;
  inputfilename2 = *++argv;
  inputfilename3 = *++argv;
  outputfilename = *++argv;
  
  inputfile3=efopen(inputfilename3,"r");
  readParmtop(inputfile3);
  fclose(inputfile3);
  numatom=AP.NATOM;
  coord=(double *)ecalloc(sizeof(double),numatom*3);
  velo =(double *)ecalloc(sizeof(double),numatom*3);  

  inputfile2=efopen(inputfilename2,"r");
  numbody=RMD_readRB_numbody(inputfile2,RBC);
  RBC.numbody =numbody;
  RBC.RB=(RigidBody *)ecalloc(sizeof(RigidBody),numbody);
  RBC.indexofABAcyc=(int *)ecalloc(sizeof(int),numbody);
  RBD_raedRBdata(inputfile2,RBC);
  fclose(inputfile2);

  inputfile1=efopen(inputfilename1,"r");
  getline(&line,&len,inputfile1);
  getline(&line,&len,inputfile1);
  io_scanconf(inputfile1,numatom,coord,'x');
  io_scanconf(inputfile1,numatom,velo, 'x');
  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      velo[i*3+j]=22.45*velo[i*3+j];
  fclose(inputfile1);

  Rotmat   = (double *)ecalloc(sizeof(double),numbody*3*3);
  Transmat = (double *)ecalloc(sizeof(double),numbody*6*6);
  spvelo   = (double *)ecalloc(sizeof(double),numbody*6);
  dtheta   = (double *)ecalloc(sizeof(double),numbody);
  crd_org  = (double *)ecalloc(sizeof(double),3);
  velo_com = (double *)ecalloc(sizeof(double),3);

  for (i=0;i<numbody;++i) {
    nparent=RBC.RB[i].num_parent_body-1;
    natom1=RBC.RB[i].num_origin-1;
    natom2=RBC.RB[i].num_origin+1-1;
    numatom_body=RBC.RB[i].numatom_body+1;
    if (i==0)
      natom3=-1;
    else
      natom3=RBC.RB[nparent].num_terminal[0]-1;
    coord_body=(double *)ecalloc(sizeof(double),numatom_body*3);
    velo_body=(double *)ecalloc(sizeof(double),numatom_body*3);
    Rmat      =(double *)ecalloc(sizeof(double),3*3);
    RBD_trans_body_coordinate(natom1,natom2,natom3,numatom_body,coord,Rmat,coord_body);
    for (j=0;j<3;++j)
      for (k=0;k<3;++k)
	Rotmat[i*9+j*3+k]=Rmat[j*3+k];
    /********************************************************************************/
    /* if (i>0) {								    */
    /*   natomi=RBC.RB[i].num_origin-1;						    */
    /*   natomj=RBC.RB[nparent].num_origin-1;					    */
    /*   Rotmati = (double *)ecalloc(sizeof(double),3*3);			    */
    /*   Rotmatj = (double *)ecalloc(sizeof(double),3*3);			    */
    /*   for (j=0;j<3;++j){							    */
    /* 	for (k=0;k<3;++k){							    */
    /* 	  Rotmati[j*3+k]=Rotmat[nparent*9+j*3+k];				    */
    /* 	  Rotmatj[j*3+k]=Rotmat[i*9+j*3+k];					    */
    /* 	}									    */
    /*   }									    */
    /*   Transmati_j = (double *)ecalloc(sizeof(double),6*6);			    */
    /*   RBD_set_trans_mat(natomi,natomj,coord,Rotmati,Rotmatj,Transmati_j);	    */
    /*   for (j=0;j<6;++j)							    */
    /* 	for (k=0;k<6;++k)							    */
    /* 	  Transmat[i*36+j*6+k]=Transmati_j[j*6+k];				    */
    /*   free(Rotmati);								    */
    /*   free(Rotmatj);								    */
    /* }									    */
    /********************************************************************************/
    mass =(double *)ecalloc(sizeof(double),numatom_body);
    for (j=0;j<numatom_body;++j) {
      mass[j]=AP.AMASS[j+natom1];
      for (k=0;k<3;++k) {
	coord_body[j*3+k] = coord[(j+natom1)*3+k];
	velo_body[j*3+k]  = velo[(j+natom1)*3+k];
      }
    }
    Inertia = (double *)ecalloc(sizeof(double),3*3);
    RBD_setcom(coord_body,mass,numatom_body,crd_org);
    RBD_setcom(velo_body,mass,numatom_body,velo_com);
    RBD_setInertia(coord_body,crd_org,mass,numatom_body,Inertia);
    //    velo_lab  = (double *)ecalloc(sizeof(double),3);
    //    velo_body = (double *)ecalloc(sizeof(double),3);
    mom       = (double *)ecalloc(sizeof(double),3);
    /**************************************/
    /* for (j=0;j<3;++j)		  */
    /*   velo_lab[j]=velo[natom1*3+j];	  */
    /* mvmult(velo_lab,Rmat,velo_body,3); */
    /**************************************/
    //    free(velo_lab);
    for (j=0;j<3;++j)
      mom[j]=0.0;
    for (j=0;j<numatom_body;++j) {
      coord_body_temp = (double *)ecalloc(sizeof(double),3);
      velo_temp       = (double *)ecalloc(sizeof(double),3);
      for (k=0;k<3;++k)
	coord_body_temp[k]=coord_body[j*3+k]-crd_org[k];
      for (k=0;k<3;++k)
	velo_temp[k]=velo_body[j*3+k]-velo_com[k]/*-velo[natom1*3+k]*/;
      //      mvmult(velo_temp,Rmat,velo_body,3);
      cv = (double *)ecalloc(sizeof(double),3);
      outprod(coord_body_temp,velo_temp,cv);
      for (k=0;k<3;++k)
	mom[k]+=mass[j]*cv[k];
      free(coord_body_temp);
      free(velo_temp);
      free(cv);
    }
    det = Inertia[0]*(Inertia[4]*Inertia[8]-Inertia[5]*Inertia[7])
         -Inertia[1]*(Inertia[3]*Inertia[8]-Inertia[6]*Inertia[5])
         +Inertia[2]*(Inertia[3]*Inertia[7]-Inertia[6]*Inertia[4]);
    InvInertia = (double *)ecalloc(sizeof(double),3*3);
    invm(Inertia,InvInertia,3);
    og = (double *)ecalloc(sizeof(double),3);
    mvmult(InvInertia,mom,og,3);
    free(Rmat);
    free(Inertia);
    free(InvInertia);
    free(mom);
    for (j=0;j<3;++j) {
      spvelo[i*6+j]=og[j];
      spvelo[i*6+j+3]=velo_body[j];
    }
    free(og);
    //    free(velo_body);
    if (i>0) {
      spvelo_pre = (double *)ecalloc(sizeof(double),6);
      theta      = (double *)ecalloc(sizeof(double),3);
      TRmat      = (double *)ecalloc(sizeof(double),3);
      TRtheta    = (double *)ecalloc(sizeof(double),3);
      /*********************************************************/
      /* Tspvelo      = (double *)ecalloc(sizeof(double),6);   */
      /* TTransmati_j = (double *)ecalloc(sizeof(double),6*6); */
      /*********************************************************/
      for (j=0;j<6;++j)
      	spvelo_pre[j]=spvelo[nparent*6+j];
      for (j=0;j<3;++j)
	theta[j] = spvelo[i*6+j]-spvelo_pre[j];
      //      mtrans(Rmat,TRmat,3);
      mvmult(Rmat,theta,TRtheta,3);
      dtheta[i]=TRtheta[2];
      free(spvelo_pre);
      free(theta);
      free(TRmat);
      free(TRtheta);
      /***********************/
      /* free(Tspvelo);	     */
      /* free(Transmati_j);  */
      /* free(TTransmati_j); */
      /***********************/
    }
    else {
      omeg_term=(double *)ecalloc(sizeof(double),3);
      velo_term=(double *)ecalloc(sizeof(double),3);
      mvmult(Rmat,spvelo,omeg_term,3);
      mvmult(Rmat,velo,velo_term,3);
    }
    free(coord_body);
    free(velo_body);
    free(mass);
  }
  
  outputfile=efopen(outputfilename,"w");
  for (i=0;i<3;++i)
    fprintf(outputfile,"%lf\n",omeg_term[i]);
  for (i=0;i<3;++i)
    fprintf(outputfile,"%lf\n",velo_term[i]);
  for (i=1;i<numbody;++i)
    fprintf(outputfile,"%lf\n",dtheta[i]);
  fclose(outputfile);

  free(RBC.RB);
  free(coord);
  free(velo);  
  free(mass);
  free(Rotmat);
  free(crd_org);
  free(velo_com);
  free(spvelo);
  free(dtheta);
  free(omeg_term);
  free(velo_term);
  
  return 0;
}
