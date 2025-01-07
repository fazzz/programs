
#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>

#include "PTL.h"
#include "FFL.h"
#include "EF.h"

#include "netcdf_mineL.h"

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}

char *NAME_ENERGY_TERMS_L[9]= {
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
    enc_def_var((nc_id_MCD->ncid),NAME_ENERGY_TERMS_L[i],NC_DOUBLE,NDIMS_ENERGY,nc_id_MCD->dimids_ene,&(nc_id_MCD->ene_term_varid[i]));
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
  if((c=nc_put_vara_double((nc_id_MCD.ncid),(nc_id_MCD.trj_varid),&nc_id_MCD.start_trj[0],&nc_id_MCD.count_trj[0],&crd_nc[0][0])))
    ERR(c);

  nc_id_MCD.start_ene[0]=numstep;
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[0]),&nc_id_MCD.start_ene[0],&ene.p_t)))
      ERR(c);
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[1]),&nc_id_MCD.start_ene[0],&ene.p_e_t)))
      ERR(c);
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[2]),&nc_id_MCD.start_ene[0],&ene.p_LJ_t)))
      ERR(c);
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[3]),&nc_id_MCD.start_ene[0],&ene.p_e_14_t)))
      ERR(c);
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[4]),&nc_id_MCD.start_ene[0],&ene.p_LJ_14_t)))
      ERR(c);
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[5]),&nc_id_MCD.start_ene[0],&ene.p_d_t)))
      ERR(c);
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[6]),&nc_id_MCD.start_ene[0],&ene.p_a_t)))
      ERR(c);
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[7]),&nc_id_MCD.start_ene[0],&ene.p_b_t)))
      ERR(c);
  if((nc_put_var1_double((nc_id_MCD.ncid),(nc_id_MCD.ene_term_varid[8]),&nc_id_MCD.start_ene[0],&drest)))
      ERR(c);

}

int myncL_put_crd_MCD(struct my_netcdf_out_id_MCD nc_id_MCD,
			 int numstep,double crd_nc[MAXATOM][3]) {
  int c;
  int i;

  nc_id_MCD.start_trj[0]=numstep;
  if((c=nc_put_vara_double((nc_id_MCD.ncid),(nc_id_MCD.trj_varid),&nc_id_MCD.start_trj[0],&nc_id_MCD.count_trj[0],&crd_nc[0][0])))
    ERR(c);

  
}


int myncL_open_inq_get_MCD(char *outfilename,int numatom,
			  int numini,int interval,int numfin,
			  struct my_netcdf_out_id_MCD *nc_id_MCD,double ***trj,double **ene){
  int i,j,k;

  enc_open(outfilename,/*NC_NOWRITE*/NC_SHARE,&(nc_id_MCD->ncid));

  enc_inq_varid(nc_id_MCD->ncid,COORDINATE,&(nc_id_MCD->trj_varid));

  for (i=0;i<NENERGY_TERMS;++i)
    enc_inq_varid(nc_id_MCD->ncid,NAME_ENERGY_TERMS_L[i],&(nc_id_MCD->ene_term_varid[i]));

  nc_id_MCD->count_trj[0] = 1;
  nc_id_MCD->count_trj[1] = numatom;
  nc_id_MCD->count_trj[2] = 2;
  nc_id_MCD->start_trj[1] = 0;
  nc_id_MCD->start_trj[2] = 0;
  nc_id_MCD->count_ene[0] = 1;

  k=0;
  for (i = numini; i < numfin; i+=interval)  {
    nc_id_MCD->start_trj[0] = i;
    nc_get_vara_double(nc_id_MCD->ncid,nc_id_MCD->trj_varid,    
		       nc_id_MCD->start_trj,nc_id_MCD->count_trj,
		       &trj[j][0][0]);                          
    nc_id_MCD->start_ene[0] = i;
    for (j=0;j<8;++j)
      nc_get_vara_double(nc_id_MCD->ncid,nc_id_MCD->ene_term_varid[j],         
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
    nc_get_vara_double(nc_id_MCD->ncid,nc_id_MCD->trj_varid,
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

  enc_inq_varid(nc_id_MCD->ncid,NAME_ENERGY_TERMS_L[index_eneterm],&(nc_id_MCD->ene_term_varid[index_eneterm]));

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

