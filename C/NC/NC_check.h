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
			      double criteria, int *ncmap
			      );

int *make_native_contact_list_aa(int *numnc, double *cord, 
				 int numatom, double criteria,
				 int *ncmap, int HMDE );

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

#endif
