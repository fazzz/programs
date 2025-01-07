#ifndef INCLUDE_PO
#define INCLUDE_PO

#include <glib.h>

#define YES 1
#define NO  0

#define ON 1
#define OFF  0

/*************************/
/* #define numatomtp 28	 */
/* #define numatomtp2 50 */
/* #define numatomtp3 10 */
/* #define numatomtp4 4	 */
/*************************/

#define numatomtp 40
#define numatomtp2 62
#define numatomtp3 17
#define numatomtp4 5

#define ON 1
#define OFF 0

//int make_bp(int numatom,int *bp,int **bp_f , int *numb, char *name_atom_list, int aromaflag);

int make_bp_v2(int numatom,int *bp,int **bp_f , int *numb, char *name_atom_list,char *name_res_list, int aromaflag, int flagtoponon);

int close_ring_v2(int numatom,int **bp_f , int *numb, char *name_atom_list, char *name_res_list);

int make_bp(int numatom,int *bp,int **bp_f , int *numb, char *name_atom_list, int aromaflag, int flagtoponon);

int order_bp(int *bp, int **bp_f, int numatom, int *numb);

void sea_p_by_name(char name[4], int nb, int *bp,int ap[2], int apflag,GHashTable *bplist,GHashTable *bplist_c2,GHashTable *bplist_c3,GHashTable *bplist_c4,GHashTable *bplist_c5,int numatom, char *name_atom_list);

int close_ring(int numatom,int **bp_f , int *numb, char *name_atom_list);

int make_nb_matrix(int **bpairs, int *numb, int depth,
		   int **matrix, int numatom);

int set_nb_pairs(int **matrix, int numatom, 
		 int **pairex1_5, int **pair1_4, 
		 int *num1_5, int *num14);

int ch_aromatic_topo(int **bp_f, int *numb,
		     int numatom,
		     char *name_atom_list);

int ch_imso_topo(int **bp_f, int *numb,
		 int numatom,
		 char *name_atom_list);

int ch_imdo_topo(int **bp_f, int *numb,
		 int numatom,
		 char *name_atom_list);

int ch_PRO_topo(int **bp_f, int *numb,
		int numatom,
		char *name_atom_list);

int ch_PRO_topo2(int **bp_f, int *numb,
		 int numatom,
		 char *name_atom_list);

int ch_ARG(int **bp_f, int *numb,
	     int numatom,
	   char *name_atom_list);

#endif

