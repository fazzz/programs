
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>

#include "netcdf_mine.h"
#include "bestfit.h"
#include "QUA.h"
#include "EF.h"
#include "PT.h"

#include "f2c.h"
#include "clapack.h"

void mmult(double *coordA, double *coordB, double matrix[3][3], int numatom,double *mass);
void fomKmat(double Kmat[4][4], double mat[3][3]);
void transCentroid(double *coordA, double *coordB, double *mass, int numatom);
void transMotion(double *velo, double *mass, int numatom);

double bestfit(double *coord_ref,double *coord_tag,double *velo_tag,double *mass,int numatom,int flag) {
  int i,j;
  double coord_bestfit[4],coord_tag_dummy[4],velo_bestfit[4],velo_tag_dummy[4];
  double mat[3][3];
  double K[4][4],Kdummy[16],rmsd=0.0;
  double q[4],eigenvalue[4],work[16];
  long int info=0,lwork=16,n=4;

  transCentroid(coord_ref,coord_tag,mass,numatom);
  mmult(coord_ref,coord_tag,mat,numatom,mass);
  fomKmat(K,mat);
 
  for (i=0;i<4;++i)
    eigenvalue[i]=0.0;
  for (i=0;i<16;++i)
    work[i]=0.0;
  
  Kdummy[0]=K[0][0];  Kdummy[4]=K[0][1];  Kdummy[8]=K[0][2];  Kdummy[12]=K[0][3];
  Kdummy[1]=K[1][0];  Kdummy[5]=K[1][1];  Kdummy[9]=K[1][2];  Kdummy[13]=K[1][3];
  Kdummy[2]=K[2][0];  Kdummy[6]=K[2][1];  Kdummy[10]=K[2][2]; Kdummy[14]=K[2][3];
  Kdummy[3]=K[3][0];  Kdummy[7]=K[3][1];  Kdummy[11]=K[3][2]; Kdummy[15]=K[3][3];

  dsyev_("V", "U", &n, Kdummy, &n, eigenvalue, work, &lwork, &info);
  for (i=0;i<4;++i)
    q[i]=Kdummy[12+i];
  for (i=0;i<numatom;++i) {
    coord_tag_dummy[0]=0.0;
    if (flag == 'b')
      velo_tag_dummy[0]=0.0;
    for (j=1;j<4;++j)  {
      coord_tag_dummy[j]=coord_tag[i*3+j-1];
      if (flag == 'b')
	velo_tag_dummy[j]=velo_tag[i*3+j-1];
    }
    qua_rot(coord_tag_dummy,q,coord_bestfit);
    if (flag == 'b')
      qua_rot(velo_tag_dummy,q,velo_bestfit);
    for (j=0;j<3;++j) {
      coord_tag[i*3+j]=coord_bestfit[j+1];
      if (flag == 'b')
	velo_tag[i*3+j]=velo_bestfit[j+1];
    }
  }

  if (flag == 'b')
    transMotion(velo_tag,mass,numatom);
  
  for (i=0;i<numatom*3;++i)
    rmsd += (coord_tag[i]-coord_ref[i])*(coord_tag[i]-coord_ref[i]);      
  rmsd = sqrt(abs(rmsd)/numatom);

  return rmsd;
}

void bf_trajectry(int numatom, int numstep, double *mass,double *coord_ref,double *rmsd_trj,char *trjname,char *velname, char *traj_bfname,char *velo_bfname,int flagv,int flaga,int flago, int flagc){
  int i,j,k,s;
  double *crd_tag,*vel_tag;
  double coord_bestfit[3],*coord_tag_old;
  double rmsd=0.0;
  double omega[3],angMom[3];
  FILE *trj,*vel,*traj_bf,*velo_bf;

  char *line;
  size_t len=0;
  static long int m=3,n=3,lda=3,info,piv[3],lwork=3;
  static double work[3];
  double Inertia[9], Inertia_dummy[3][3];

  trj=efopen(trjname,"r");
  if (flagv=='b')
    vel=efopen(velname,"r");
  if (flago=='o') {
    traj_bf=efopen(traj_bfname,"w");
    if (flagv=='b')
      velo_bf=efopen(velo_bfname,"w");
  }
  
  crd_tag=(double *)gcemalloc(sizeof(double)*numatom*3);
  vel_tag=(double *)gcemalloc(sizeof(double)*numatom*3);

  if (flaga=='a')
    getline(&line,&len,trj);
  if (flaga=='a' && flagv =='b')
    getline(&line,&len,vel);
  for (s=0;s<numstep;++s) {
    if (flagc=='c') {
      io_scancaconf(trj,crd_tag);
      if (flagv=='b')
	io_scancaconf(vel,vel_tag);
    }
    else {
      io_scanconf(trj,numatom,crd_tag,'x');
      if (flagv=='b')
	io_scanconf(vel,numatom,vel_tag,'x');
    }
    rmsd_trj[s]=bestfit(coord_ref,crd_tag,vel_tag,mass,numatom,flagv);
    
    if (flagv=='b') {
      for (i=0;i<3;++i)
	angMom[i] = 0.0;
      for (i=0;i<numatom;++i){
	angMom[0]+= mass[i]*(crd_tag[i*3+2]*vel_tag[i*3+1]-crd_tag[i*3+1]*vel_tag[i*3+2]);
	angMom[1]+= mass[i]*(crd_tag[i*3+0]*vel_tag[i*3+2]-crd_tag[i*3+2]*vel_tag[i*3+0]);
	angMom[2]+= mass[i]*(crd_tag[i*3+1]*vel_tag[i*3+0]-crd_tag[i*3+0]*vel_tag[i*3+1]);
      }

      for (i=0;i<9;++i)
	Inertia[i]=0.0;
      for (i=0;i<numatom;++i){
	Inertia[0]+=mass[i]*crd_tag[i*3+1]*crd_tag[i*3+1]+mass[i]*crd_tag[i*3+2]*crd_tag[i*3+2];
	Inertia[1]-=mass[i]*crd_tag[i*3+0]*crd_tag[i*3+1];
	Inertia[2]-=mass[i]*crd_tag[i*3+0]*crd_tag[i*3+2];

	Inertia[3]-=mass[i]*crd_tag[i*3+1]*crd_tag[i*3+0];
	Inertia[4]+=mass[i]*crd_tag[i*3+0]*crd_tag[i*3+0]+mass[i]*crd_tag[i*3+2]*crd_tag[i*3+2];
	Inertia[5]-=mass[i]*crd_tag[i*3+1]*crd_tag[i*3+2];
	
	Inertia[6]-=mass[i]*crd_tag[i*3+2]*crd_tag[i*3+0];
	Inertia[7]-=mass[i]*crd_tag[i*3+2]*crd_tag[i*3+1];
	Inertia[8]+=mass[i]*crd_tag[i*3+0]*crd_tag[i*3+0]+mass[i]*crd_tag[i*3+1]*crd_tag[i*3+1];
      }

      dgetrf_(&m,&n,Inertia,&lda,piv,&info);
      dgetri_(&n,Inertia,&lda,piv,work,&lwork,&info);

      k=0;
      for (i=0;i<3;++i){
	for (j=0;j<3;++j){
	  Inertia_dummy[i][j]=Inertia[k];
	  ++k;
	}
      }

      for (i=0;i<3;++i)
	omega[i]=0.0;
      for (i=0;i<3;++i){
	for (j=0;j<3;++j){
	  omega[i]+=Inertia_dummy[i][j]*angMom[j];
	}
      }

      for (i=0;i<numatom;++i)  {
	vel_tag[i*3+0] -= crd_tag[i*3+1]*omega[2]-crd_tag[i*3+2]*omega[1];
	vel_tag[i*3+1] -= crd_tag[i*3+2]*omega[0]-crd_tag[i*3+0]*omega[2];
	vel_tag[i*3+2] -= crd_tag[i*3+0]*omega[1]-crd_tag[i*3+1]*omega[0];
      }

      for (i=0;i<3;++i)
	angMom[i] = 0.0;
      for (i=0;i<numatom;++i){
	angMom[0]+= mass[i]*(crd_tag[i*3+2]*vel_tag[i*3+1]-crd_tag[i*3+1]*vel_tag[i*3+2]);
	angMom[1]+= mass[i]*(crd_tag[i*3+0]*vel_tag[i*3+2]-crd_tag[i*3+2]*vel_tag[i*3+0]);
	angMom[2]+= mass[i]*(crd_tag[i*3+1]*vel_tag[i*3+0]-crd_tag[i*3+0]*vel_tag[i*3+1]);
      }
    }
    
    if (flago=='o') {
      for (j=0;j<numatom;++j) {
	for (k=0;k<3;++k) 
	  fprintf(traj_bf,"%12.8lf ",crd_tag[j*3+k]);
	fprintf(traj_bf,"\n");
      }
      fprintf(traj_bf,"\n");
      
      if (flagv=='b') {
	for (j=0;j<numatom;++j) {
	  for (k=0;k<3;++k) 
	    fprintf(velo_bf,"%12.8lf ",vel_tag[j*3+k]);
	  fprintf(velo_bf,"\n");
	}
	fprintf(velo_bf,"\n");
      }
    }
  }

  fclose(trj);
  if (flagv=='b')
    fclose(vel);
  if (flago=='o') {
    fclose(traj_bf);
    if (flagv=='b')
      fclose(velo_bf);
  }
}

void mmult(double *coordA, double *coordB, double matrix[3][3], int numatom,double *mass) {
  int i,j,k;

  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      matrix[i][j] = 0.0;
  for (i=0;i<3;++i)
    for (j=0;j<3;++j)
      for (k=0;k<numatom;++k)
	matrix[i][j] += coordA[k*3+i]*coordB[k*3+j]*mass[k];
}

void fomKmat(double Kmat[4][4], double mat[3][3]) {
  Kmat[0][0] =  mat[0][0]+mat[1][1]+mat[2][2];
  Kmat[1][1] =  mat[0][0]-mat[1][1]-mat[2][2];
  Kmat[2][2] = -mat[0][0]+mat[1][1]-mat[2][2];
  Kmat[3][3] = -mat[0][0]-mat[1][1]+mat[2][2];

  Kmat[0][1] = mat[1][2]-mat[2][1];
  Kmat[0][2] = mat[2][0]-mat[0][2];
  Kmat[0][3] = mat[0][1]-mat[1][0];

  Kmat[1][2] = mat[0][1]+mat[1][0];
  Kmat[2][3] = mat[1][2]+mat[2][1];
  Kmat[1][3] = mat[2][0]+mat[0][2];

  Kmat[1][0] = Kmat[0][1];
  Kmat[2][0] = Kmat[0][2];
  Kmat[2][1] = Kmat[1][2];

  Kmat[3][0] = Kmat[0][3];
  Kmat[3][1] = Kmat[1][3];
  Kmat[3][2] = Kmat[2][3];

}

void transCentroid(double *coordA, double *coordB, double *mass, int numatom) {
  int i,j;
  double cm_A[3],cm_B[3];
  double summass;

  summass=0.0;

  for (j=0;j<3;++j)  {
    cm_A[j] = 0.0;
    cm_B[j] = 0.0; 
  }

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j) {
      cm_A[j] += coordA[i*3+j]*mass[i];
      cm_B[j] += coordB[i*3+j]*mass[i]; 
    }
    summass += mass[i];
  }

 for (i=0;i<3;++i)  {
   cm_A[i] = cm_A[i]/(summass);
   cm_B[i] = cm_B[i]/(summass); 
 }

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j)  {
      coordA[i*3+j] -= cm_A[j]; 
      coordB[i*3+j] -= cm_B[j]; 
    }
  }

}

void transMotion(double *velo, double *mass, int numatom) {
  int i,j;
  double v_COM[3];
  double summass;

  summass=0.0;
  for (j=0;j<3;++j) 
    v_COM[j] = 0.0;

  for (i=0;i<numatom;++i) {
    for (j=0;j<3;++j)
      v_COM[j] += velo[i*3+j]*mass[i];
    summass += mass[i];
  }
  
 for (i=0;i<3;++i)
   v_COM[i] = v_COM[i]/(summass);

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j)
      velo[i*3+j] -= v_COM[j]; 

}

void ave_coord(double *coord_ref , FILE *inputfile, int numatom, int numstep, int flagc) {
  int i,j,s;
  double *crd_tag;
  
  crd_tag=(double *)gcemalloc(sizeof(double)*numatom*3);
  for (s=0;s<numstep;++s) {
    if (flagc=='c')
      io_scancaconf(inputfile,crd_tag);
    else
      io_scanconf(inputfile,numatom,crd_tag,'x');
    for (i=0;i<numatom;++i)
      for (j=0;j<3;++j)
	coord_ref[i*3+j]=(s*coord_ref[i*3+j]+crd_tag[i*3+j])/(s+1);
  }
}

void bf_trajectry_nc(int numatom, int numstep, 
		     double *mass,double *coord_ref,
		     double *rmsd_trj,char *trjname,char *traj_bfname,
		     int MODE,int flago, int IOMODE){
  int i,j,k,s;
  int numatomp;
  double *crd_tag,*crd_ref,*massp;
  double rmsd=0.0;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD,nc_id_MCD_bf;
  struct my_netcdf_out_id_AMBER nc_id_MD,nc_id_MD_bf;

  if (IOMODE==AMBER) {
    numstep=mync_get_present_step_AMBER(trjname,&nc_id_MD);
    mync_create_def_AMBER(traj_bfname,numatom,&nc_id_MD_bf);
  }
  else {
    numstep=mync_get_present_step_MCD(trjname,&nc_id_MCD);
    mync_create_def_MCD(traj_bfname,numatom,&nc_id_MCD_bf);
  }

  crd_tag=(double *)gcemalloc(sizeof(double)/**numatom*3*/);
  crd_ref=(double *)gcemalloc(sizeof(double)/**numatom*3*/);
  massp=(double *)gcemalloc(sizeof(double));
  for (s=0;s<numstep;++s) {
    if (IOMODE==AMBER)
      mync_open_inq_get_sh_AMBER(trjname,numatom,s,1,s+1,&nc_id_MD,crd_nc);
    else
      mync_open_inq_get_sh_MCD(trjname,numatom,s,1,s+1,&nc_id_MCD,crd_nc);
    k=0;
    for (i=0;i<numatom;++i) {
      if (MODE==AA) {
	crd_tag=(double *)gcerealloc(crd_tag,sizeof(double)*(k+1)*3);
	crd_ref=(double *)gcerealloc(crd_ref,sizeof(double)*(k+1)*3);
	massp=(double *)gcerealloc(massp,sizeof(double)*(k+1));
	for (j=0;j<3;++j) crd_tag[k*3+j]=crd_nc[i][j];
	for (j=0;j<3;++j) crd_ref[k*3+j]=coord_ref[i*3+j];
	massp[k]=mass[i];
	++k;
      }
      else if (MODE==CA) {
	if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
	  crd_tag=(double *)gcerealloc(crd_tag,sizeof(double)*(k+1)*3);
	  crd_ref=(double *)gcerealloc(crd_ref,sizeof(double)*(k+1)*3);
	  massp=(double *)gcerealloc(massp,sizeof(double)*(k+1));
	  for (j=0;j<3;++j) crd_tag[k*3+j]=crd_nc[i][j];
	  for (j=0;j<3;++j) crd_ref[k*3+j]=coord_ref[i*3+j];
	  massp[k]=mass[i];
	  ++k;
	}
      }
      else if (MODE==HV) {
	if (strncmp(AP.IGRAPH[i],"H",1)!=0) {
	  crd_tag=(double *)gcerealloc(crd_tag,sizeof(double)*(k+1)*3);
	  crd_ref=(double *)gcerealloc(crd_ref,sizeof(double)*(k+1)*3);
	  massp=(double *)gcerealloc(massp,sizeof(double)*(k+1));
	  for (j=0;j<3;++j) crd_tag[k*3+j]=crd_nc[i][j];
	  for (j=0;j<3;++j) crd_ref[k*3+j]=coord_ref[i*3+j];
	  massp[k]=mass[i];
	  ++k;
	}
      }
    }
    numatomp=k;
    rmsd_trj[s]=bestfit_nc(crd_ref,crd_tag,massp,numatomp);
    k=0;
    for (i=0;i<numatom;++i) {
      if (MODE==AA) {
	for (j=0;j<3;++j)
	  crd_nc[i][j]=crd_tag[k*3+j];
	++k;
      }
      else if (MODE==CA) {
	if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
	  for (j=0;j<3;++j)
	    crd_nc[i][j]=crd_tag[k*3+j];
	  ++k;
	}
	else {
	  for (j=0;j<3;++j)
	    crd_nc[i][j]=0.0;
	}
      }
      else if (IOMODE==HV) {
	if (strncmp(AP.IGRAPH[i],"H",1)==0) {
	  for (j=0;j<3;++j)
	    crd_nc[i][j]=crd_tag[k*3+j];
	  ++k;
	}
	else {
	  for (j=0;j<3;++j)
	    crd_nc[i][j]=0.0;
	}
      }
    }
    
    if (IOMODE==AMBER)
      mync_put_crd_AMBER(nc_id_MD_bf,s,crd_nc);
    else 
      mync_put_crd_MCD(nc_id_MCD_bf,s,crd_nc);
    /*********************************/
    /* if (IOMODE==AMBER)	     */
    /*   nc_close((nc_id_MD.ncid));  */
    /* else 			     */
    /*   nc_close((nc_id_MCD.ncid)); */
    /*********************************/
  }
  /************************************/
  /* if (IOMODE==AMBER)		      */
  /*   nc_close((nc_id_MD_bf.ncid));  */
  /* else 			      */
  /*   nc_close((nc_id_MCD_bf.ncid)); */
  /************************************/
}

void bf_trajectry_ncb(int numatom, int numstep, int interval,
		      double *mass,double *coord_ref,
		      double *rmsd_trj,char *trjname,char *traj_bfname,
		      int MODE,int flago, int IOMODE){
  int i,j,k,s,ss=0;
  int numatomp;
  double *crd_tag,*crd_ref,*massp;
  double rmsd=0.0;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD,nc_id_MCD_bf;
  struct my_netcdf_out_id_AMBER nc_id_MD,nc_id_MD_bf;

  if (IOMODE==AMBER) {
    numstep=mync_get_present_step_AMBER(trjname,&nc_id_MD);
    mync_create_def_AMBER(traj_bfname,numatom,&nc_id_MD_bf);
  }
  else {
    numstep=mync_get_present_step_MCD(trjname,&nc_id_MCD);
    mync_create_def_MCD(traj_bfname,numatom,&nc_id_MCD_bf);
  }

  crd_tag=(double *)gcemalloc(sizeof(double)/**numatom*3*/);
  crd_ref=(double *)gcemalloc(sizeof(double)/**numatom*3*/);
  massp=(double *)gcemalloc(sizeof(double));
  for (s=0;s<numstep;++s) {
    if (IOMODE==AMBER)
      mync_open_inq_get_sh_AMBER(trjname,numatom,s,1,s+1,&nc_id_MD,crd_nc);
    else
      mync_open_inq_get_sh_MCD(trjname,numatom,s,1,s+1,&nc_id_MCD,crd_nc);

    k=0;
    for (i=0;i<numatom;++i) {
      if (MODE==AA) {
	crd_tag=(double *)gcerealloc(crd_tag,sizeof(double)*(k+1)*3);
	crd_ref=(double *)gcerealloc(crd_ref,sizeof(double)*(k+1)*3);
	massp=(double *)gcerealloc(massp,sizeof(double)*(k+1));
	for (j=0;j<3;++j) crd_tag[k*3+j]=crd_nc[i][j];
	for (j=0;j<3;++j) crd_ref[k*3+j]=coord_ref[i*3+j];
	massp[k]=mass[i];
	++k;
      }
      else if (MODE==CA) {
	if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
	  crd_tag=(double *)gcerealloc(crd_tag,sizeof(double)*(k+1)*3);
	  crd_ref=(double *)gcerealloc(crd_ref,sizeof(double)*(k+1)*3);
	  massp=(double *)gcerealloc(massp,sizeof(double)*(k+1));
	  for (j=0;j<3;++j) crd_tag[k*3+j]=crd_nc[i][j];
	  for (j=0;j<3;++j) crd_ref[k*3+j]=coord_ref[i*3+j];
	  massp[k]=mass[i];
	  ++k;
	}
      }
      else if (MODE==HV) {
	if (strncmp(AP.IGRAPH[i],"H",1)!=0) {
	  crd_tag=(double *)gcerealloc(crd_tag,sizeof(double)*(k+1)*3);
	  crd_ref=(double *)gcerealloc(crd_ref,sizeof(double)*(k+1)*3);
	  massp=(double *)gcerealloc(massp,sizeof(double)*(k+1));
	  for (j=0;j<3;++j) crd_tag[k*3+j]=crd_nc[i][j];
	  for (j=0;j<3;++j) crd_ref[k*3+j]=coord_ref[i*3+j];
	  massp[k]=mass[i];
	  ++k;
	}
      }
    }
    numatomp=k;
    rmsd_trj[s]=bestfit_nc(crd_ref,crd_tag,massp,numatomp);
    k=0;
    for (i=0;i<numatom;++i) {
      if (MODE==AA) {
	for (j=0;j<3;++j)
	  crd_nc[i][j]=crd_tag[k*3+j];
	++k;
      }
      else if (MODE==CA) {
	if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
	  for (j=0;j<3;++j)
	    crd_nc[i][j]=crd_tag[k*3+j];
	  ++k;
	}
	else {
	  for (j=0;j<3;++j)
	    crd_nc[i][j]=0.0;
	}
      }
      else if (IOMODE==HV) {
	if (strncmp(AP.IGRAPH[i],"H",1)==0) {
	  for (j=0;j<3;++j)
	    crd_nc[i][j]=crd_tag[k*3+j];
	  ++k;
	}
	else {
	  for (j=0;j<3;++j)
	    crd_nc[i][j]=0.0;
	}
      }
    }
    
    if (s%interval==0) {
      if (IOMODE==AMBER)
	mync_put_crd_AMBER(nc_id_MD_bf,ss,crd_nc);
      else 
	mync_put_crd_MCD(nc_id_MCD_bf,ss,crd_nc);
      ++ss;
    }
    /*********************************/
    /* if (IOMODE==AMBER)	     */
    /*   nc_close((nc_id_MD.ncid));  */
    /* else 			     */
    /*   nc_close((nc_id_MCD.ncid)); */
    /*********************************/
  }
  /************************************/
  /* if (IOMODE==AMBER)		      */
  /*   nc_close((nc_id_MD_bf.ncid));  */
  /* else 			      */
  /*   nc_close((nc_id_MCD_bf.ncid)); */
  /************************************/

}


void ave_coord_nc(double *coord_ref ,char *trjname,int numatom,int numstep,int MODE, int IOMODE) {
  int i,j,k,s;
  double *crd_tag;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;
  
  crd_tag=(double *)gcemalloc(sizeof(double)*numatom*3);
  for (s=0;s<numstep;++s) {
    if (IOMODE==AMBER)
      mync_open_inq_get_sh_AMBER(trjname,numatom,s,1,s+1,&nc_id_MD,crd_nc);
    else
      mync_open_inq_get_sh_MCD(trjname,numatom,s,1,s+1,&nc_id_MCD,crd_nc);
    k=0;
    for (i=0;i<numatom;++i) {
      if (MODE==AA) {
	for (j=0;j<3;++j)
	  crd_tag[k*3+j]=crd_nc[i][j];
	++k;
      }
      else if (MODE==CA) {
	if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
	  for (j=0;j<3;++j)
	    crd_tag[i*3+j]=crd_nc[i][j];
	}
	else {
	  for (j=0;j<3;++j)
	    crd_tag[i*3+j]=0.0;
	}
      }
      else if (MODE==HV) {
	if (strncmp(AP.IGRAPH[i],"H",1)!=0) {
	  for (j=0;j<3;++j)
	    crd_tag[i*3+j]=crd_nc[i][j];
	}
	else {
	  for (j=0;j<3;++j)
	    crd_tag[i*3+j]=0.0;
	}
      }
    }
    for (i=0;i<numatom;++i)
      for (j=0;j<3;++j)
	coord_ref[i*3+j]=(s*coord_ref[i*3+j]+crd_tag[i*3+j])/(s+1);
  }
}

void ave_coord_ncb(double *coord_ref ,char *trjname,int numatom,int numstep,int interval,int MODE, int IOMODE) {
  int i,j,k,s,ss=0;
  double *crd_tag;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;
  
  crd_tag=(double *)gcemalloc(sizeof(double)*numatom*3);
  for (s=0;s<numstep;++s) {
    if (IOMODE==AMBER)
      mync_open_inq_get_sh_AMBER(trjname,numatom,s,1,s+1,&nc_id_MD,crd_nc);
    else
      mync_open_inq_get_sh_MCD(trjname,numatom,s,1,s+1,&nc_id_MCD,crd_nc);
    if (s%interval==0) {
      k=0;
      for (i=0;i<numatom;++i) {
	if (MODE==AA) {
	  for (j=0;j<3;++j)
	    crd_tag[k*3+j]=crd_nc[i][j];
	  ++k;
	}
	else if (MODE==CA) {
	  if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
	    for (j=0;j<3;++j)
	      crd_tag[i*3+j]=crd_nc[i][j];
	  }
	  else {
	    for (j=0;j<3;++j)
	      crd_tag[i*3+j]=0.0;
	  }
	}
	else if (MODE==HV) {
	  if (strncmp(AP.IGRAPH[i],"H",1)!=0) {
	    for (j=0;j<3;++j)
	      crd_tag[i*3+j]=crd_nc[i][j];
	  }
	  else {
	    for (j=0;j<3;++j)
	      crd_tag[i*3+j]=0.0;
	  }
	}
      }
      for (i=0;i<numatom;++i)
	for (j=0;j<3;++j)
	  coord_ref[i*3+j]=(ss*coord_ref[i*3+j]+crd_tag[i*3+j])/(ss+1);
      ++ss;
    }
  }
}


double bestfit_nc(double *coord_ref,double *coord_tag,double *mass,int numatom) {
  int i,j;
  double coord_bestfit[4],coord_tag_dummy[4];
  double mat[3][3];
  double K[4][4],Kdummy[16],rmsd=0.0;
  double q[4],eigenvalue[4],work[16];
  long int info=0,lwork=16,n=4;

  transCentroid(coord_ref,coord_tag,mass,numatom);
  mmult(coord_ref,coord_tag,mat,numatom,mass);
  fomKmat(K,mat);
 
  for (i=0;i<4;++i)
    eigenvalue[i]=0.0;
  for (i=0;i<16;++i)
    work[i]=0.0;
  
  Kdummy[0]=K[0][0];  Kdummy[4]=K[0][1];  Kdummy[8]=K[0][2];  Kdummy[12]=K[0][3];
  Kdummy[1]=K[1][0];  Kdummy[5]=K[1][1];  Kdummy[9]=K[1][2];  Kdummy[13]=K[1][3];
  Kdummy[2]=K[2][0];  Kdummy[6]=K[2][1];  Kdummy[10]=K[2][2]; Kdummy[14]=K[2][3];
  Kdummy[3]=K[3][0];  Kdummy[7]=K[3][1];  Kdummy[11]=K[3][2]; Kdummy[15]=K[3][3];

  dsyev_("V", "U", &n, Kdummy, &n, eigenvalue, work, &lwork, &info);
  for (i=0;i<4;++i)
    q[i]=Kdummy[12+i];
  for (i=0;i<numatom;++i) {
    coord_tag_dummy[0]=0.0;
    for (j=1;j<4;++j)  {
      coord_tag_dummy[j]=coord_tag[i*3+j-1];
    }
    qua_rot(coord_tag_dummy,q,coord_bestfit);
    for (j=0;j<3;++j) {
      coord_tag[i*3+j]=coord_bestfit[j+1];
    }
  }

  for (i=0;i<numatom*3;++i)
    rmsd += (coord_tag[i]-coord_ref[i])*(coord_tag[i]-coord_ref[i]);      
  rmsd = sqrt(abs(rmsd)/numatom);

  return rmsd;
}
