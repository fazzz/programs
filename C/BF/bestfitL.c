
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <netcdf.h>

#include "netcdf_mineL.h"
#include "bestfit.h"
#include "bestfitL.h"
#include "PTL.h"

void bf_trajectry_ncbL(int numatom, int numstep, int interval,
		       double *mass,double *coord_ref,
		       double *rmsd_trj,char *trjname,char *traj_bfname,
		       int MODE,int flago, int IOMODE, int INMODE){
  int i,j,k,s,ss=0;
  int numatomp;
  double *crd_tag,*crd_ref,*massp;
  double rmsd=0.0;

  double crd_nc[MAXATOM][3];

  struct my_netcdf_out_id_MCD nc_id_MCD,nc_id_MCD_bf;
  struct my_netcdf_out_id_AMBER nc_id_MD,nc_id_MD_bf;

  if (IOMODE==AMBER) {
    numstep=myncL_get_present_step_AMBER(trjname,&nc_id_MD);
    myncL_create_def_AMBER(traj_bfname,numatom,&nc_id_MD_bf);
  }
  else {
    numstep=myncL_get_present_step_MCD(trjname,&nc_id_MCD);
    myncL_create_def_MCD(traj_bfname,numatom,&nc_id_MCD_bf);
  }

  crd_tag=(double *)gcemalloc(sizeof(double)/**numatom*3*/);
  crd_ref=(double *)gcemalloc(sizeof(double)/**numatom*3*/);
  massp=(double *)gcemalloc(sizeof(double));
  for (s=0;s<numstep;++s) {
    if (IOMODE==AMBER) myncL_open_inq_get_sh_AMBER(trjname,numatom,s,1,s+1,&nc_id_MD,crd_nc);
    else myncL_open_inq_get_sh_MCD(trjname,numatom,s,1,s+1,&nc_id_MCD,crd_nc);

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
      else if (MODE==HV) {
	if (strncmp(AP.IGRAPH[i],"H",1)!=0) {
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
      if (IOMODE==AMBER) myncL_put_crd_AMBER(nc_id_MD_bf,ss,crd_nc);
      else myncL_put_crd_MCD(nc_id_MCD_bf,ss,crd_nc);
      ++ss;
    }
  }

}

void ave_coord_ncbL(double *coord_ref ,char *trjname,int numatom,int numstep,int interval,int MODE, int IOMODE) {
  int i,j,k,s,ss=0;
  double *crd_tag;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;
  
  crd_tag=(double *)gcemalloc(sizeof(double)*numatom*3);
  for (s=0;s<numstep;++s) {
    if (IOMODE==AMBER) myncL_open_inq_get_sh_AMBER(trjname,numatom,s,1,s+1,&nc_id_MD,crd_nc);
    else myncL_open_inq_get_sh_MCD(trjname,numatom,s,1,s+1,&nc_id_MCD,crd_nc);
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
