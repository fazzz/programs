#ifndef INCLUDE_NC
#define INCLUDE_NC

#define YES 1
#define NO 0

#define ON 1
#define OFF 0

#define NC 1
#define NN 0

#define INC 0
#define EXC 1

#define criteria_NC 6.5

int *make_native_contact_list(int *numnc, double *cord, 
			      int numatom, int numres,
			      double criteria, int **ncmap
			      );

int *make_native_contact_list_aa(int *numnc, double *cord, 
				 int numatom, double criteria,
				 int **ncmap, int HMDE );

int *make_native_contact_list_aa_2(int *numncaa, int *numncres,double *cord, 
				   int numatom, int numres, double criteria,
				   int **ncmap, int HMODE);

int *make_native_contact_list_aa_3(int *numncaa, int *numncres,double *cord, 
				   int numatom, int numres, double criteria,
				   int **ncmap,  int **ncmapres,int HMODE);

int *make_native_contact_list_aa_3_nadjacent(int *numncaa, int *numncres,double *cord, 
					     int numatom, int numres, double criteria,
					     int **ncmap,  int **ncmapres,int HMODE);

int *make_native_contact_list_aa_3_nadjacent_2(int *numncaa, int *numncres,double *cord, 
					       int numatom, int numres, double criteria,
					       int **ncmap,  int **ncmapres,int HMODE);

double count_native_contact(int numnc, double *cord, 
			    int numatom, int numres,
			    int *indexncb, double criteria, int HMDE);

/********************************************************************/
/* int count_native_contact(int numnc, double *cord, 		    */
/* 			 int numatom, int numres,		    */
/* 			 int *indexncb, double criteria);	    */
/********************************************************************/

//int trace_fraction_q_value();

int check_within_3_neibor_res(int num1, int num2);

double count_native_contact_aa(int numnc, double *crd, 
			       int numatom,int numres,
			       int **ncmap, double *cradii_natatt);

double count_native_contact_hybrid(int numncaa, double *crd, 
				   int numatom,int numres,
				   int **ncmap, double *cradii_natatt);

double count_native_contact_ca(int numnc,double *crd, 
			       int numres,
			       int **ncmapres,double *cradii_natatt);

int *make_native_contact_list_aa_wnibnum(int *numncaa, int *numncres,double *cord, 
					 int numatom, int numres, double criteria,
					 int **ncmap,  int **ncmapres,int HMODE, int nibnum);

#endif
