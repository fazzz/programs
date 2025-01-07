
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>

#include "PTL.h"
#include "FFL.h"
#include "EF.h"
#include "SBFF.h"

#include "netcdf_mine.h"
//#include "netcdf_mineL.h"

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}

char *NAME_ENERGY_TERMS[9]= {
  "energy_total",
  "energy_ele_sta",
  "energy_vdW",
  "energy_ele_sta_1_4",
  "energy_vdW_1_4",
  "energy_torsion",
  "energy_angle",
  "energy_bond",
  "energy_rest"
};

int myncL_create_def_MCD(char *outfilename,int numatom,
			struct my_netcdf_out_id_MCD *nc_id_MCD){
  int i;


  enc_create(outfilename,NC_SHARE,&(nc_id_MCD->ncid));

  enc_def_dim(nc_id_MCD->ncid,XYZ,NXYZ,&(nc_id_MCD->xyz_dimid));
  enc_def_dim(nc_id_MCD->ncid,MOLECULE,numatom,&(nc_id_MCD->molecule_dimid));
  enc_def_dim(nc_id_MCD->ncid,REC_NAME,NC_UNLIMITED,&(nc_id_MCD->rec_dimid));

  nc_id_MCD->dimids_trj[0] = (nc_id_MCD->rec_dimid);
  nc_id_MCD->dimids_trj[1] = (nc_id_MCD->molecule_dimid);
  nc_id_MCD->dimids_trj[2] = (nc_id_MCD->xyz_dimid);

  nc_id_MCD->dimids_ene[0] = (nc_id_MCD->rec_dimid);

  enc_def_var((nc_id_MCD->ncid),COORDINATE,NC_DOUBLE,NDIMS_TRJ,nc_id_MCD->dimids_trj,&(nc_id_MCD->trj_varid));
  nc_put_att_text((nc_id_MCD->ncid),(nc_id_MCD->trj_varid),UNITS,strlen(TRJ_UNIT),TRJ_UNIT);

  for (i=0;i<NENERGY_TERMS;++i) {
    enc_def_var((nc_id_MCD->ncid),NAME_ENERGY_TERMS[i],NC_DOUBLE,NDIMS_ENERGY,nc_id_MCD->dimids_ene,&(nc_id_MCD->ene_term_varid[i]));
    nc_put_att_text((nc_id_MCD->ncid),(nc_id_MCD->ene_term_varid[i]),UNITS,strlen(ENE_UNIT),ENE_UNIT);
  }

  nc_enddef((nc_id_MCD->ncid));

  nc_id_MCD->count_trj[0] = 1;
  nc_id_MCD->count_trj[1] = numatom;
  nc_id_MCD->count_trj[2] = 3;
  nc_id_MCD->start_trj[1] = 0;
  nc_id_MCD->start_trj[2] = 0;
  nc_id_MCD->count_ene[0] = 1;

}

int myncL_put_crd_ene_MCD(struct my_netcdf_out_id_MCD nc_id_MCD,
			 int numstep,double crd_nc[MAXATOM][3],
			 struct potential ene, double drest) {
  int c;
  int i;

  nc_id_MCD.start_trj[0]=numstep;
  if((c=nc_put_vara_double((nc_id_MCD.ncid),(nc_id_MCD.trj_varid),&nc_id_MCD.start_trj,&nc_id_MCD.count_trj,&crd_nc[0][0])))
    ERR(c);

  nc_id_MCD.start_ene[0]=numstep;
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[0]),&nc_id_MCD.start_ene,&ene.p_t)))
      ERR(c);
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[1]),&nc_id_MCD.start_ene,&ene.p_e_t)))
      ERR(c);
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[2]),&nc_id_MCD.start_ene,&ene.p_LJ_t)))
      ERR(c);
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[3]),&nc_id_MCD.start_ene,&ene.p_e_14_t)))
      ERR(c);
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[4]),&nc_id_MCD.start_ene,&ene.p_LJ_14_t)))
      ERR(c);
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[5]),&nc_id_MCD.start_ene,&ene.p_d_t)))
      ERR(c);
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[6]),&nc_id_MCD.start_ene,&ene.p_a_t)))
      ERR(c);
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[7]),&nc_id_MCD.start_ene,&ene.p_b_t)))
      ERR(c);
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[8]),&nc_id_MCD.start_ene,&drest)))
      ERR(c);

}

int myncL_put_crd_MCD(struct my_netcdf_out_id_MCD nc_id_MCD,
			 int numstep,double crd_nc[MAXATOM][3]) {
  int c;
  int i;

  nc_id_MCD.start_trj[0]=numstep;
  if((c=nc_put_vara_double((nc_id_MCD.ncid),(nc_id_MCD.trj_varid),&nc_id_MCD.start_trj,&nc_id_MCD.count_trj,&crd_nc[0][0])))
    ERR(c);

  
}


int myncL_open_inq_get_MCD(char *outfilename,int numatom,
			  int numini,int interval,int numfin,
			  struct my_netcdf_out_id_MCD *nc_id_MCD,double ***trj,double **ene){
  int i,j,k;

  enc_open(outfilename,/*NC_NOWRITE*/NC_SHARE,&(nc_id_MCD->ncid));

  enc_inq_varid(nc_id_MCD->ncid,COORDINATE,&(nc_id_MCD->trj_varid));

  for (i=0;i<NENERGY_TERMS;++i)
    enc_inq_varid(nc_id_MCD->ncid,NAME_ENERGY_TERMS[i],&(nc_id_MCD->ene_term_varid[i]));

  nc_id_MCD->count_trj[0] = 1;
  nc_id_MCD->count_trj[1] = numatom;
  nc_id_MCD->count_trj[2] = 2;
  nc_id_MCD->start_trj[1] = 0;
  nc_id_MCD->start_trj[2] = 0;
  nc_id_MCD->count_ene[0] = 1;

  k=0;
  for (i = numini; i < numfin; i+=interval)  {
    nc_id_MCD->start_trj[0] = i;
    nc_get_vara_double(nc_id_MCD->ncid,COORDINATE,
		       nc_id_MCD->start_trj,nc_id_MCD->count_trj,
		       &trj[j][0][0]);
    nc_id_MCD->start_ene[0] = i;
    for (j=0;j<8;++j)
      nc_get_vara_double(nc_id_MCD->ncid,NAME_ENERGY_TERMS[j],
			 nc_id_MCD->start_ene,nc_id_MCD->count_ene,&ene[k][j]);
    ++k;
  }

  encclose(nc_id_MCD->ncid);

}

int myncL_open_inq_get_trj_MCD(char *outfilename,int numatom,
			  int numini,int interval,int numfin,
			  struct my_netcdf_out_id_MCD *nc_id_MCD,double ***trj){
  int i,j,k;

  enc_open(outfilename,NC_NOWRITE,&(nc_id_MCD->ncid));

  enc_inq_varid(nc_id_MCD->ncid,COORDINATE,&(nc_id_MCD->trj_varid));

  nc_id_MCD->count_trj[0] = 1;
  nc_id_MCD->count_trj[1] = numatom;
  nc_id_MCD->count_trj[2] = 2;
  nc_id_MCD->start_trj[1] = 0;
  nc_id_MCD->start_trj[2] = 0;
  nc_id_MCD->count_ene[0] = 1;

  k=0;
  for (i = numini; i < numfin; i+=interval)  {
    nc_id_MCD->start_trj[0] = i;
    nc_get_vara_double(nc_id_MCD->ncid,COORDINATE,
		       nc_id_MCD->start_trj,nc_id_MCD->count_trj,
		       &trj[j][0][0]);
    nc_id_MCD->start_ene[0] = i;
    ++k;
  }

  encclose(nc_id_MCD->ncid);
}

int myncL_open_inq_get_sh_MCD(char *outfilename,int numatom,
			     int numini,int interval,int numfin,
			     struct my_netcdf_out_id_MCD *nc_id_MCD,
			     double crd_nc[MAXATOM][3]){
  int i,j,k;
  int c;

  enc_open(outfilename,NC_NOWRITE,&(nc_id_MCD->ncid));

  enc_inq_varid(nc_id_MCD->ncid,COORDINATE,&(nc_id_MCD->trj_varid));

  nc_id_MCD->count_trj[0] = 1;
  nc_id_MCD->count_trj[1] = numatom;
  nc_id_MCD->count_trj[2] = 3;
  nc_id_MCD->start_trj[1] = 0;
  nc_id_MCD->start_trj[2] = 0;

  k=0;
  for (i = numini; i < numfin; i+=interval)  {
    nc_id_MCD->start_trj[0] = i;
    //    nc_id_MCD->count_trj[0] = i+1;
    if (c=nc_get_vara_double(nc_id_MCD->ncid,nc_id_MCD->trj_varid,nc_id_MCD->start_trj,nc_id_MCD->count_trj,&crd_nc[0][0])) {
      ERR(c);
    }
    ++k;
  }

  encclose(nc_id_MCD->ncid);

}

int myncL_open_inq_get_ene_MCD(char *infilename,
			      int numini,int interval,int numfin,int index_eneterm,
			      struct my_netcdf_out_id_MCD *nc_id_MCD,
			      double *ene){
  int i,j,k;
  int c;

  enc_open(infilename,NC_NOWRITE,&(nc_id_MCD->ncid));

  enc_inq_varid(nc_id_MCD->ncid,NAME_ENERGY_TERMS[index_eneterm],&(nc_id_MCD->ene_term_varid[index_eneterm]));

  nc_id_MCD->count_ene[0] = 1;

  k=0;
  for (i=numini;i<numfin;i+=interval)  {
    nc_id_MCD->start_ene[0] = i;
    if (c=nc_get_vara_double(nc_id_MCD->ncid,
			     nc_id_MCD->ene_term_varid[index_eneterm],
			     nc_id_MCD->start_ene,nc_id_MCD->count_ene,ene))
	ERR(c);
    ++k;
  }

  encclose(nc_id_MCD->ncid);

}

int myncL_get_present_step_MCD(char *infilename,
			      struct my_netcdf_out_id_MCD *nc_id_MCD){
  int elength;
  int c,unlimdimid;
  size_t length;

  enc_open(infilename,NC_NOWRITE,&(nc_id_MCD->ncid));
  if(c=nc_inq_unlimdim(nc_id_MCD->ncid,&unlimdimid))
    ERR(c);

  if((c=nc_inq_dimlen(nc_id_MCD->ncid,unlimdimid,&length)))
    ERR(c);
  elength=length;

  encclose(nc_id_MCD->ncid);
  
  return elength;
}

int myncL_get_numatom_MCD(char *infilename,
			 struct my_netcdf_out_id_MCD *nc_id_MCD){
  int numatom;
  int c,molecule_dimid;
  size_t length;

  enc_open(infilename,NC_NOWRITE,&(nc_id_MCD->ncid));
  if(c=nc_inq_dimid(nc_id_MCD->ncid,MOLECULE,&molecule_dimid))
    ERR(c);

  if((c=nc_inq_dimlen(nc_id_MCD->ncid,molecule_dimid,&length)))
    ERR(c);
  numatom=length;

  encclose(nc_id_MCD->ncid);
  
  return numatom;
}

///////////////////////////////////////////////////////////////////////////////////////////

char *NAME_SBAA_ENERGY_TERMS[6]= {
  "energy_total",
  "energy_cnb_t",
  "energy_nnb_t",
  "energy_d_t",
  "energy_a_t",
  "energy_b_t",
};

int myncL_create_def_SBAAMCD(char *outfilename,int numatom,
			    struct my_netcdf_out_id_SBAAMCD *nc_id_MCD){
  int i;

  enc_create(outfilename,NC_SHARE,&(nc_id_MCD->ncid));

  enc_def_dim(nc_id_MCD->ncid,XYZ,NXYZ,&(nc_id_MCD->xyz_dimid));
  enc_def_dim(nc_id_MCD->ncid,MOLECULE,numatom,&(nc_id_MCD->molecule_dimid));
  enc_def_dim(nc_id_MCD->ncid,REC_NAME,NC_UNLIMITED,&(nc_id_MCD->rec_dimid));

  nc_id_MCD->dimids_trj[0] = (nc_id_MCD->rec_dimid);
  nc_id_MCD->dimids_trj[1] = (nc_id_MCD->molecule_dimid);
  nc_id_MCD->dimids_trj[2] = (nc_id_MCD->xyz_dimid);

  nc_id_MCD->dimids_ene[0] = (nc_id_MCD->rec_dimid);

  enc_def_var((nc_id_MCD->ncid),COORDINATE,NC_DOUBLE,NDIMS_TRJ,nc_id_MCD->dimids_trj,&(nc_id_MCD->trj_varid));
  nc_put_att_text((nc_id_MCD->ncid),(nc_id_MCD->trj_varid),UNITS,strlen(TRJ_UNIT),TRJ_UNIT);

  for (i=0;i<NENERGY_SBAA_TERMS;++i) {
    enc_def_var((nc_id_MCD->ncid),NAME_SBAA_ENERGY_TERMS[i],NC_DOUBLE,NDIMS_ENERGY,nc_id_MCD->dimids_ene,&(nc_id_MCD->ene_term_varid[i]));
    nc_put_att_text((nc_id_MCD->ncid),(nc_id_MCD->ene_term_varid[i]),UNITS,strlen(ENE_UNIT),ENE_UNIT);
  }

  nc_enddef((nc_id_MCD->ncid));

  nc_id_MCD->count_trj[0] = 1;
  nc_id_MCD->count_trj[1] = numatom;
  nc_id_MCD->count_trj[2] = 3;
  nc_id_MCD->start_trj[1] = 0;
  nc_id_MCD->start_trj[2] = 0;
  nc_id_MCD->count_ene[0] = 1;

}

int myncL_put_crd_ene_SBAAMCD(struct my_netcdf_out_id_SBAAMCD nc_id_MCD,
			     int numstep,double crd_nc[MAXATOM][3],
			     /*struct potential_SBAA ene*/ double p_t) {
  int c;
  int i;

  nc_id_MCD.start_trj[0]=numstep;
  if((c=nc_put_vara_double((nc_id_MCD.ncid),(nc_id_MCD.trj_varid),&nc_id_MCD.start_trj,&nc_id_MCD.count_trj,&crd_nc[0][0])))
    ERR(c);

  nc_id_MCD.start_ene[0]=numstep;
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[0]),&nc_id_MCD.start_ene,&/*ene.*/p_t)))
      ERR(c);
  /**************************************************************************************************************/
  /* if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[1]),&nc_id_MCD.start_ene,&ene.p_cnb_t))) */
  /*     ERR(c);											        */
  /* if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[2]),&nc_id_MCD.start_ene,&ene.p_nnb_t))) */
  /*     ERR(c);											        */
  /* if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[3]),&nc_id_MCD.start_ene,&ene.p_d_t)))   */
  /*     ERR(c);											        */
  /* if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[4]),&nc_id_MCD.start_ene,&ene.p_a_t)))   */
  /*     ERR(c);											        */
  /* if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[5]),&nc_id_MCD.start_ene,&ene.p_b_t)))   */
  /*     ERR(c);											        */
  /**************************************************************************************************************/

}

/************************************************************************************************************************************************************************/
/* int myncL_put_crd_ene_each_SBAAMCD(struct my_netcdf_out_id_SBAAMCD nc_id_MCD,int numstep,double crd_nc[MAXATOM][3],int index_eneterm,struct potential_SBAA ene p_t) { */
/*   int c;																			        */
/*   int i;																			        */
/* 																				        */
/*   nc_id_MCD.start_trj[0]=numstep;																        */
/*   if((c=nc_put_vara_double((nc_id_MCD.ncid),(nc_id_MCD.trj_varid),&nc_id_MCD.start_trj,&nc_id_MCD.count_trj,&crd_nc[0][0])))					        */
/*     ERR(c);																			        */
/* 																				        */
/*   nc_id_MCD.start_ene[0]=numstep;																        */
/*   if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[index_eneterm]),&nc_id_MCD.start_ene,&/\*ene.*\/p_t)))					        */
/*       ERR(c);																		        */
/* }																				        */
/************************************************************************************************************************************************************************/


int myncL_open_inq_get_ene_SBAAMCD(char *infilename,
				  int numini,int interval,int numfin,
				  struct my_netcdf_out_id_SBAAMCD *nc_id_MCD,
				  double *ene){
  int i,j,k;
  int c;

  enc_open(infilename,NC_NOWRITE,&(nc_id_MCD->ncid));

  enc_inq_varid(nc_id_MCD->ncid,NAME_ENERGY_TERMS[0],&(nc_id_MCD->ene_term_varid[0]));

  nc_id_MCD->count_ene[0] = 1;

  k=0;
  for (i=numini;i<numfin;i+=interval)  {
    nc_id_MCD->start_ene[0] = i;
    if (c=nc_get_vara_double(nc_id_MCD->ncid,
			     nc_id_MCD->ene_term_varid[0],
			     nc_id_MCD->start_ene,nc_id_MCD->count_ene,ene))
	ERR(c);
    ++k;
  }

  encclose(nc_id_MCD->ncid);

}

int myncL_get_present_step_SBAAMCD(char *infilename,
				  struct my_netcdf_out_id_SBAAMCD *nc_id_MCD){
  int elength;
  int c,unlimdimid;
  size_t length;

  enc_open(infilename,NC_NOWRITE,&(nc_id_MCD->ncid));
  if(c=nc_inq_unlimdim(nc_id_MCD->ncid,&unlimdimid))
    ERR(c);

  if((c=nc_inq_dimlen(nc_id_MCD->ncid,unlimdimid,&length)))
    ERR(c);
  elength=length;

  encclose(nc_id_MCD->ncid);
  
  return elength;
}


///////////////////////////////////////////////////////////////////////////////////////////

int myncL_create_def_AMBER(char *outfilename,int numatom,
			  struct my_netcdf_out_id_AMBER *nc_id_MD){

  enc_create(outfilename,NC_SHARE,&(nc_id_MD->ncid));

  enc_def_dim(nc_id_MD->ncid,"spatial",3,&(nc_id_MD->spatial_dimid));
  enc_def_dim(nc_id_MD->ncid,"atom",numatom,&(nc_id_MD->atom_dimid));
  enc_def_dim(nc_id_MD->ncid,"frame",NC_UNLIMITED,&(nc_id_MD->rec_dimid));

  nc_id_MD->dimids_trj[0] = (nc_id_MD->rec_dimid);
  nc_id_MD->dimids_trj[1] = (nc_id_MD->atom_dimid);
  nc_id_MD->dimids_trj[2] = (nc_id_MD->spatial_dimid);

  enc_def_var((nc_id_MD->ncid),"coordinates",NC_DOUBLE,3,nc_id_MD->dimids_trj,&(nc_id_MD->trj_varid));
  nc_put_att_text((nc_id_MD->ncid),(nc_id_MD->trj_varid),"units",strlen("angstrom"),"angstrom");

  nc_enddef((nc_id_MD->ncid));

  nc_id_MD->count_trj[0] = 1;
  nc_id_MD->count_trj[1] = numatom;
  nc_id_MD->count_trj[2] = 3;
  nc_id_MD->start_trj[1] = 0;
  nc_id_MD->start_trj[2] = 0;

}

int myncL_get_present_step_AMBER(char *infilename,
				struct my_netcdf_out_id_AMBER *nc_id_MD){
  int elength;
  int c,unlimdimid;
  size_t length;

  enc_open(infilename,NC_NOWRITE,&(nc_id_MD->ncid));
  if(c=nc_inq_unlimdim(nc_id_MD->ncid,&unlimdimid))
    ERR(c);

  if((c=nc_inq_dimlen(nc_id_MD->ncid,unlimdimid,&length)))
    ERR(c);
  elength=length;

  encclose(nc_id_MD->ncid);
  
  return elength;
}

int myncL_open_inq_get_sh_AMBER(char *outfilename,int numatom,
			       int numini,int interval,int numfin,
			       struct my_netcdf_out_id_AMBER *nc_id_MD,
			       double crd_nc[MAXATOM][3]){
  int i,j,k;
  int c;

  enc_open(outfilename,NC_NOWRITE,&(nc_id_MD->ncid));

  enc_inq_varid(nc_id_MD->ncid,"coordinates",&(nc_id_MD->trj_varid));

  nc_id_MD->count_trj[0] = 1;
  nc_id_MD->count_trj[1] = numatom;
  nc_id_MD->count_trj[2] = 3;
  nc_id_MD->start_trj[1] = 0;
  nc_id_MD->start_trj[2] = 0;

  k=0;
  for (i = numini; i < numfin; i+=interval)  {
    nc_id_MD->start_trj[0] = i;

    if (c=nc_get_vara_double(nc_id_MD->ncid,nc_id_MD->trj_varid,nc_id_MD->start_trj,nc_id_MD->count_trj,&crd_nc[0][0])) {
      ERR(c);
    }
    ++k;
  }

  encclose(nc_id_MD->ncid);

}

int myncL_put_crd_AMBER(struct my_netcdf_out_id_AMBER nc_id_MD,
		       int numstep,double crd_nc[MAXATOM][3]) {
  int c;
  int i;

  nc_id_MD.start_trj[0]=numstep;
  if((c=nc_put_vara_double((nc_id_MD.ncid),(nc_id_MD.trj_varid),&nc_id_MD.start_trj,&nc_id_MD.count_trj,&crd_nc[0][0])))
    ERR(c);

  
}

double *myncL_get_trj_aw(char *inputfilename,int MODE,int IOMODE,
			int numatom,int *numatomp,int *numstep){
  int s,i,j,k,l;

  double *crd,*mass;
  double *traj;

  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;
  crd=(double *)gcemalloc(sizeof(double)/**numatom*3*/);
  mass=(double *)gcemalloc(sizeof(double));
  traj=(double *)gcemalloc(sizeof(double));

  if (IOMODE==AMBER) *numstep=myncL_get_present_step_AMBER(inputfilename,&nc_id_MD);
  else *numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MCD);

  for (i=0;i<*numstep;++i) {
    if (IOMODE==AMBER)
      myncL_open_inq_get_sh_AMBER(inputfilename,numatom,i,1,i+1,&nc_id_MD,crd_nc);
    else
      myncL_open_inq_get_sh_MCD(inputfilename,numatom,i,1,i+1,&nc_id_MCD,crd_nc);

    l=0;
    for (j=0;j<numatom;++j) {
      if (MODE==AA) {
	crd=(double *)gcerealloc(crd,sizeof(double)*(l+1)*3);	
	mass=(double *)gcerealloc(mass,sizeof(double)*(l+1));
	for (k=0;k<3;++k) crd[l*3+k]=crd_nc[j][k];
	mass[l]=AP.AMASS[j];
	++l;
      }
      else if (MODE==CA) {
	if (strncmp(AP.IGRAPH[j],"CA",2)==0) {
	  crd=(double *)gcerealloc(crd,sizeof(double)*(l+1)*3);
	  mass=(double *)gcerealloc(mass,sizeof(double)*(l+1));
	  for (k=0;k<3;++k) crd[l*3+k]=crd_nc[j][k];
	  mass[l]=AP.AMASS[j];
	  ++l;
	}
      }
      else if (MODE==HV) {
	if (strncmp(AP.IGRAPH[j],"H",1)!=0) {
	  crd=(double *)gcerealloc(crd,sizeof(double)*(l+1)*3);
	  mass=(double *)gcerealloc(mass,sizeof(double)*(l+1));
	  for (k=0;k<3;++k) crd[l*3+k]=crd_nc[j][k];
	  mass[l]=AP.AMASS[j];
	  ++l;
	}
      }
    }
    *numatomp=l;
    traj=(double *)gcerealloc(traj,sizeof(double)*(*numstep)*(*numatomp)*3);
    for (j=0;j<*numatomp;++j)  {
      for (k=0;k<3;++k) {
	traj[i*(*numatomp)*3+j*3+k] = crd[j*3+k]*sqrt(mass[j]);
      }
    }
  }

  return traj;
}

double *myncL_get_trj_aw_b(char *inputfilename,int inMODE,int outMODE,int IOMODE,
			   int numatom,int *numatomp,int *numstep){
  int s,i,j,k,l,m;

  int numatompin;
  double *crd,*mass,*traj;
  double crd_nc[MAXATOM][3];
  struct my_netcdf_out_id_MCD nc_id_MCD;
  struct my_netcdf_out_id_AMBER nc_id_MD;

  mass=(double *)gcemalloc(sizeof(double));

  if ((inMODE==CA && outMODE==AA) || (inMODE==CA && outMODE==HV) 
      || (inMODE==HV && outMODE==AA) || (inMODE==HV && outMODE==CA) ) 
    printf("error\n");

  if (IOMODE==AMBER) *numstep=myncL_get_present_step_AMBER(inputfilename,&nc_id_MD);
  else *numstep=myncL_get_present_step_MCD(inputfilename,&nc_id_MCD);

  l=0;
  m=0;
  for (i=0;i<numatom;++i) {
    if (outMODE==AA) {
      mass=(double *)gcerealloc(mass,sizeof(double)*(l+1));
      mass[l]=AP.AMASS[i];
      ++l;
    }
    else if (outMODE==CA) {
      if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
	mass=(double *)gcerealloc(mass,sizeof(double)*(l+1));
	mass[l]=AP.AMASS[i];
	++l;
      }
    }
    else if (outMODE==HV){
      if (strncmp(AP.IGRAPH[i],"H",1)!=0) {
	mass=(double *)gcerealloc(mass,sizeof(double)*(l+1));
	mass[l]=AP.AMASS[i];
	++l;
      }
    }
    if (inMODE==AA) ++m;
    else if (inMODE==CA) if (strncmp(AP.IGRAPH[i],"CA",2)==0) ++m;
    else if (inMODE==HV) if (strncmp(AP.IGRAPH[i],"H",1)!=0) ++m;
  }
  *numatomp=l;
  numatompin=m;
  crd=(double *)gcemalloc(sizeof(double)*(*numatomp)*3);
  traj=(double *)gcemalloc(sizeof(double)*(*numstep)*(*numatomp)*3);

  for (i=0;i<*numstep;++i) {
    if (IOMODE==AMBER)
      myncL_open_inq_get_sh_AMBER(inputfilename,numatompin,i,1,i+1,&nc_id_MD,crd_nc);
    else
      myncL_open_inq_get_sh_MCD(inputfilename,numatompin,i,1,i+1,&nc_id_MCD,crd_nc);

    l=0;
    for (j=0;j<numatompin;++j) {
      if (inMODE==AA) {
	if (outMODE==AA) {
	  for (k=0;k<3;++k) crd[j*3+k]=crd_nc[j][k];
	}
	if (outMODE==CA) {
	  if (strncmp(AP.IGRAPH[i],"CA",2)==0) {
	    for (k=0;k<3;++k) crd[l*3+k]=crd_nc[j][k];
	    ++l;
	  }
	}
	if (outMODE==HV) {
	  if (strncmp(AP.IGRAPH[i],"H",1)!=0) {
	    for (k=0;k<3;++k) crd[l*3+k]=crd_nc[j][k];
	    ++l;
	  }
	}
      }
      else
	for (k=0;k<3;++k) crd[j*3+k]=crd_nc[j][k];
    }      
      
    for (j=0;j<*numatomp;++j) 
      for (k=0;k<3;++k)
	traj[i*(*numatomp)*3+j*3+k] = crd[j*3+k]*sqrt(mass[j]);
  }
    
  return traj;
}

///////////////////////////////////////////////////////////////////////////////////////////
