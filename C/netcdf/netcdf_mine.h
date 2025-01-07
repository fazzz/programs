#ifndef INCLUDE_netcdf
#define INCLUDE_netcdf

#include "FF.h"

#define NDIMS_TRJ 3
#define NDIMS_ENERGY 1

#define NXYZ 3
#define NENERGY_TERMS 9

// dimensions
#define XYZ "xyz"
#define MOLECULE "molecule"

#define REC_NAME "step"

#define UNITS "units"
#define TRJ_UNIT "angstroum"
#define ENE_UNIT "kcal_mol"

// variables
#define COORDINATE "coordinate"

//#define MAXATOM 300
//#define MAXATOM 3000
#define MAXATOM 10000

#define MD 0
#define AMBER 1

#define AA 0
#define CA 1
#define HV 2

struct my_netcdf_out_id_MCD {
  int ncid;
  int xyz_dimid,molecule_dimid,rec_dimid;
  int trj_varid;

  int ene_term_varid[NENERGY_TERMS];

  size_t start_trj[NDIMS_TRJ],count_trj[NDIMS_TRJ];
  size_t start_ene[1],count_ene[1];

  int dimids_trj[NDIMS_TRJ],dimids_ene[NDIMS_ENERGY];
};

int mync_create_def_MCD(char *outfilename,int numatom,
			struct my_netcdf_out_id_MCD *nc_id_MCD);

int mync_put_crd_ene_MCD(struct my_netcdf_out_id_MCD nc_id_MCD,
			 int numstep,
			 /*double *crd_nc*/
			 double crd_nc[MAXATOM][3],
			 struct potential ene, double drest);

int mync_put_crd_MCD(struct my_netcdf_out_id_MCD nc_id_MCD,
		     int numstep,
		     double crd_nc[MAXATOM][3]);

int mync_open_inq_get_MCD(char *outfilename,int numatom,
			  int numini,int interval,int numfin,
			  struct my_netcdf_out_id_MCD *nc_id_MCD,
			  double ***trj,double **ene);

int mync_open_inq_get_trj_MCD(char *outfilename,int numatom,
			      int numini,int interval,int numfin,
			      struct my_netcdf_out_id_MCD *nc_id_MCD,
			      double ***trj);

int mync_open_inq_get_sh_MCD(char *outfilename,int numatom,
			     int numini,int interval,int numfin,
			     struct my_netcdf_out_id_MCD *nc_id_MCD,
			     double crd_nc[MAXATOM][3]);

int mync_open_inq_get_ene_MCD(char *infilename,
			      int numini,int interval,int numfin,int index_eneterm,
			      struct my_netcdf_out_id_MCD *nc_id_MCD,
			      double *ene);

int mync_get_present_step_MCD(char *infilename,
			      struct my_netcdf_out_id_MCD *nc_id_MCD);

int mync_get_numatom_MCD(char *infilename,
			 struct my_netcdf_out_id_MCD *nc_id_MCD);

///////////////////////////////////////////////////////////////////////////////////////

#define NENERGY_SBAA_TERMS 1

struct my_netcdf_out_id_SBAAMCD {
  int ncid;
  int xyz_dimid,molecule_dimid,rec_dimid;
  int trj_varid;

  int ene_term_varid[NENERGY_SBAA_TERMS];

  size_t start_trj[NDIMS_TRJ],count_trj[NDIMS_TRJ];
  size_t start_ene[1],count_ene[1];

  int dimids_trj[NDIMS_TRJ],dimids_ene[NDIMS_ENERGY];
};

int mync_create_def_SBAAMCD(char *outfilename,int numatom,
			     struct my_netcdf_out_id_SBAAMCD *nc_id_MCD);

int mync_put_crd_ene_SBAAMCD(struct my_netcdf_out_id_SBAAMCD nc_id_MCD,
			     int numstep,double crd_nc[MAXATOM][3],
			     /*struct potential_SBAA ene*/double p_t);

int mync_open_inq_get_ene_SBAAMCD(char *infilename,
				  int numini,int interval,int numfin,
				  struct my_netcdf_out_id_SBAAMCD *nc_id_MCD,double *ene);

int mync_get_present_step_SBAAMCD(char *infilename,
				  struct my_netcdf_out_id_SBAAMCD *nc_id_MCD);

////////////////////////////////////////////////////////////////////

struct my_netcdf_out_id_AMBER {
  int ncid;
  int spatial_dimid,atom_dimid,rec_dimid;
  int trj_varid;

  size_t start_trj[NDIMS_TRJ],count_trj[NDIMS_TRJ];

  int dimids_trj[NDIMS_TRJ];
};

int mync_create_def_AMBER(char *outfilename,int numatom,
			  struct my_netcdf_out_id_AMBER *nc_id_MD);

int mync_get_present_step_AMBER(char *infilename,
				struct my_netcdf_out_id_AMBER *nc_id_MD);

int mync_open_inq_get_sh_AMBER(char *outfilename,int numatom,
			       int numini,int interval,int numfin,
			       struct my_netcdf_out_id_AMBER *nc_id_MD,
			       double crd_nc[MAXATOM][3]);

int mync_put_crd_AMBER(struct my_netcdf_out_id_AMBER nc_id_MD,
		       int numstep,double crd_nc[MAXATOM][3]);

double *mync_get_trj_aw(char *inputfilename,int MODE,int IOMODE,
			int numatom,int *numatomp,int *numstep);

#endif
