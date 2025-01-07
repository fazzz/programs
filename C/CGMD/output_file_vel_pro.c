#include <stdio.h>
#include <math.h>
#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "physics.h"
#include "MD.h"

void output_file_vel_pro(FILE *output_v) {
	int i,alpha;
	int j,k,a=0;
	double Energy_kinetic_db;

        /****************************************************************/
        /* calc_velo2();					        */
	/* 							        */
	/* Energy_kinetic_db = 0.0;				        */
	/* for (i=0;i<prot.DOF;++i){				        */
	/*   for (j=0;j<clust[i].num_atom_clust;++j){		        */
	/*     for (alpha=0;alpha<3;++alpha){			        */
	/*       Energy_kinetic_db += 0.5*clust[i].mass_clust[j]        */
	/* 	*prot.velo[a][alpha]				        */
	/* 	*prot.velo[a][alpha];				        */
	/*     }						        */
	/*     ++a;						        */
	/*   }							        */
	/* }							        */
	/* Energy_kinetic_db = Energy_kinetic_db/(4.18407*100);	        */
        /****************************************************************/
	
	for(i=0; i<prot.num_atom; ++i){
	  for(alpha=0; alpha<3; ++alpha){
	    fprintf(output_v,"%12.8lf  " ,prot.velo[i][alpha]);
	  }
	  fprintf(output_v,"\n");
	}

	fprintf(output_v,"\n\n");
}

void output_file_dihed_vel_pro(FILE *output_v) {
  int i;

  for(i=0; i<prot.DOF; ++i){
    fprintf(output_v,"%12.8lf  \n" ,clust[i].ddihedang[0]);
  }
  
  fprintf(output_v,"\n\n");
}
