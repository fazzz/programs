#ifndef INCLUDE_NC_count
#define INCLUDE_NC_count

#define NCratio_default 1.2
#define NCext_default 1.0

double NCcount_native_contact_wratio(int numnc, double *crd, 
				     int numatom,int numres,
				     int **ncmap,int **ncmap_res,
				     double *cradii_natatt, double nc_ratio);

double NCcount_native_contact_wext(int numnc, double *crd, 
				   int numatom,int numres,
				   int **ncmap,int **ncmap_res,
				   double *cradii_natatt, double nc_ext);

double NCcount_native_contact_AA_wratio(int numnc, double *crd, 
					int numatom,int numres,
					int **ncmap,int **ncmap_res,
					double *cradii_natatt, double nc_ratio);

double NCcount_native_contact_AA_wext(int numnc, double *crd, 
				      int numatom,int numres,
				      int **ncmap, int **ncmap_res,
				      double *cradii_natatt, double nc_ext);

#endif


