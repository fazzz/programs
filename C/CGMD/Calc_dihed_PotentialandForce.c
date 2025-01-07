#include <stdio.h>
#include <math.h>
#include "gener.h"
#include "ABA.h" // #include "ABA_multimer.h"
#include "MD.h"
#include "force.h"
#include "ParmTop.h"

double pick_dihed_one_clust2(int nNumAom1,
	                      int nNumAom2,
	                      int nNumAom3,
	                      int nNumAom4,
			     int nNumDihed
			     // double sn,
			     //double cs
			     );

double pick_dihed_one_clust4(int nNumAom1,
	                      int nNumAom2,
	                      int nNumAom3,
	                      int nNumAom4
			     );

//void Calc_restraintForce(double dihedang,int nNumDihed,int nNumClut );

int which_rest(int dihed[5]);

// 2 面角相互作用の計算を行う関数
void Calc_dihed_Potential_Force(double *eig)
{
	int i;
	int nNumDihed;
	int nNumClut;
	int nNumDihedType;

	int flag;

	double P_dihed[3000];
	double N_dihed[3000];

	double V_N_now;
	double n_dihed_now;
	double isou_dihed_now;
	double dihedang;

	double P_dihed_debug;
	double P_dihed_debug2=0.0;
	FILE *debug,*data;

	potential_pro.p_dihedt = 0.0;

        /*********************************/
        /* debug=fopen("debug.txt","a"); */
        /*********************************/

	for (nNumClut=0;nNumClut<prot.DOF;++nNumClut)
	{
		P_dihed[nNumClut] = 0.0;
		N_dihed[nNumClut] = 0.0;
	}

        for (nNumClut=0;nNumClut<prot.DOF;++nNumClut){
	  potential_pro.p_dihedc[nNumClut] = 0.0;
	  clust[nNumClut].f_c.f_dihed = 0.0;
	}

	if (MASSIVEOUT==ON && (nNumStep%out_put_steps_thomo)==0){
	  data=fopen("dihed.txt","a");
	  fprintf(data,"%d  ",nNumStep);
	}

	for (nNumDihed=0;nNumDihed<prot.nNumDihedALL/*100*//*390*//*500*//*800*//*380*//*AP.NPHIH/*30*/;++nNumDihed) {
	  flag=ON;
          for (i=0;i<inpnum;++i) {
	    if (nNumDihed == inpindex[i]) {
	      if (MASSIVEOUT==ON && (nNumStep%out_put_steps_thomo)==0){
	  	fprintf(data,"0.0  ");
	      }
	      flag=OFF;
	      break;
	    }
	  }
	  if (flag==ON) {
	  nNumDihedType = atom_dihed_pair[nNumDihed/*][*/*6+4]-1;
	  nNumClut = atom_dihed_pair[nNumDihed/*][*/*6+5]-1/*2*//*1110*//*1112*/;

	    //	  if (nNumClut != 68) {
	    dihedang = pick_dihed_one_clust4(atom_dihed_pair[nNumDihed/*][*/*6+0]-1,
					     atom_dihed_pair[nNumDihed/*][*/*6+1]-1,
					     atom_dihed_pair[nNumDihed/*][*/*6+2]-1,
					     atom_dihed_pair[nNumDihed/*][*/*6+3]-1);

            /**********************************************************************************************************************************/
            /* fprintf(debug,"%d %d %e %e %e  \n",nNumDihed,nNumDihedType,AP.PK[nNumDihedType],AP.PN[nNumDihedType],AP.PHASE[nNumDihedType]); */
            /**********************************************************************************************************************************/
	    /******************************/
	    /* printf("%lf \n",dihedang); */
	    /******************************/
	  
	    n_dihed_now = AP.PN[nNumDihedType];
	    V_N_now =AP.PK[nNumDihedType];
	    isou_dihed_now = AP.PHASE[nNumDihedType];

	    //  printf("yes-i1\n");
	    
	    //	    debug=fopen("debug_force.txt","a");

	    if (PEPCAACCMODE == 0 ) {
	    potential_pro.p_dihedc[nNumClut] += AP.PK[nNumDihedType]*(1.0+cos(AP.PN[nNumDihedType]*dihedang-AP.PHASE[nNumDihedType]));
	    clust[nNumClut].f_c.f_dihed += -AP.PK[nNumDihedType]*4.18407*100.0
	      *(sin(AP.PN[nNumDihedType]*dihedang-AP.PHASE[nNumDihedType])*AP.PN[nNumDihedType]);

	    //fprintf(debug,"%d %d %e %e %e \n",nNumStep,nNumDihed,AP.PK[nNumDihedType]*(1.0+cos(AP.PN[nNumDihedType]*dihedang-AP.PHASE[nNumDihedType])),-AP.PK[nNumDihedType]*4.18407*100.0*(sin(AP.PN[nNumDihedType]*dihedang-AP.PHASE[nNumDihedType])*AP.PN[nNumDihedType]),dihedang);

	    }
	    else if (PEPCAACCMODE == 1 ) {
	    potential_pro.p_dihedc[nNumClut] += AP.PK[nNumDihedType]*(1.0+cos(AP.PN[nNumDihedType]*dihedang-AP.PHASE[nNumDihedType]));
	    clust[nNumClut].f_c.f_dihed += -AP.PK[nNumDihedType]*4.18407*100.0*(sin(AP.PN[nNumDihedType]*dihedang-AP.PHASE[nNumDihedType])*AP.PN[nNumDihedType]);
	    }
	    else {	      
	      //  printf("yes-i1\n");
	      potential_pro.p_dihedc[nNumClut] = AP.PK[nNumDihedType]*(1.0+cos(AP.PN[nNumDihedType]*dihedang-AP.PHASE[nNumDihedType]))*eig[nNumDihed]*fact;
	    clust[nNumClut].f_c.f_dihed = -AP.PK[nNumDihedType]*4.18407*100.0
	      *(sin(AP.PN[nNumDihedType]*dihedang-AP.PHASE[nNumDihedType])*AP.PN[nNumDihedType])*eig[nNumDihed]*fact;
	    }

	    if (MASSIVEOUT==ON && (nNumStep%out_put_steps_thomo)==0){
	      fprintf(data,"%e  ",AP.PK[nNumDihedType]*(1.0+cos(AP.PN[nNumDihedType]*dihedang-AP.PHASE[nNumDihedType])));
	    }

	    //fclose(debug);

	    //  printf("yes-i2\n");

	    /**************************************************************************************************************************************************/
            /* printf("%d \n t=%d  n=%d V=%lf isou=%lf \n",nNumDihed,nNumDihedType,AP.PN[nNumDihedType],AP.PK[nNumDihedType],AP.PHASE[nNumDihedType]);	      */
	    /* printf("%d ",nNumClut);															      */
	    /* printf("ang=%lf ",dihedang);														      */
	    /* printf("pot=%lf ",AP.PK[nNumDihedType]*(1.0+cos(AP.PN[nNumDihedType]*dihedang-AP.PHASE[nNumDihedType])));				      */
	    /* printf("foc=%lf\n",-AP.PK[nNumDihedType]*4.18407*100.0											      */
	    /*   *(sin(AP.PN[nNumDihedType]*dihedang-AP.PHASE[nNumDihedType])*AP.PN[nNumDihedType]));							      */
            /**************************************************************************************************************************************************/


	  /***********************************************************************/
	  /* potential_pro.p_dihedc[nNumClut] += cos(dihedang);		   */
	  /* clust[nNumClut].f_c.f_dihed += -4.18407*100.0*sin(dihedang);	   */
	  /***********************************************************************/
	  //	  }
	  /*****************************************************************************************************************/
	  /* P_dihed_debug = AP.PK[nNumDihedType]*(1.0+cos(AP.PN[nNumDihedType]*dihedang-AP.PHASE[nNumDihedType]));	     */
	  /*****************************************************************************************************************/
	  /**********************************************************************************/
	  /* fprintf(debug,"%d %12.10e %12.10e\n",nNumDihed,P_dihed_debug,dihedang);	      */
	  /* P_dihed_debug2 += P_dihed_debug;						      */
	  /**********************************************************************************/
	  }
	  }
        /******************************************************************/
        /* for (nNumClut=0;nNumClut<prot.DOF;++nNumClut){		  */
	/*   potential_pro.p_dihedc[nNumClut] = P_dihed[nNumClut];	  */
	/*   clust[nNumClut].f_c.f_dihed = N_dihed[nNumClut];		  */
	/* }								  */
        /******************************************************************/

	/**********************************************************/
        /************************************************************************/
        /* fprintf(debug,"\n  ");					        */
	/* /\* 							  *\/	        */
        /************************************************************************/
	/*************************/
        /******************/
        /* fclose(debug); */
        /******************/
        /*************************/
        /**********************************************************/
	if (MASSIVEOUT==ON && (nNumStep%out_put_steps_thomo)==0){
	  fprintf(data,"\n  ");
	  fclose(data);
	}
	//  printf("yes-if\n");

}

// calculation of restraint forc
void Calc_restraintforce(void/*double dihedang,int nnumdihed,int nnumclut*/) {
  int i;
  int nNumClut;
  int nNumDihed;
  double dihedang;
  //  file *db;

  potential_pro.p_restt = 0.0;

  for (nNumClut=0;nNumClut<prot.DOF;++nNumClut)	{
    potential_pro.p_rest[nNumClut] = 0.0;
    clust[nNumClut].f_c.f_rest = 0.0;
  }

  for (nNumDihed=0;nNumDihed<prot.nNumDihed_rest;++nNumDihed) {
    nNumClut = dihed_rest[nNumDihed][4]-1;

    dihedang = pick_dihed_one_clust4(dihed_rest[nNumDihed][0]-1,
				     dihed_rest[nNumDihed][1]-1,
				     dihed_rest[nNumDihed][2]-1,
				     dihed_rest[nNumDihed][3]-1);


    potential_pro.p_rest[nNumClut] = V_rest[nNumDihed]*0.5*(dihedang-theta_ref[nNumDihed])*(dihedang-theta_ref[nNumDihed]);
    //    potential_pro.p_restt += potential_pro.p_rest[nNumClut];
    clust[nNumClut].f_c.f_rest = V_rest[nNumDihed]*(dihedang-theta_ref[nNumDihed])*4.18407*100.0;
  }
}

int which_rest(int dihed[5]) 
{
  int i;
  int flag;
  
  flag = 0;
  for (i=0;i<prot.nNumDihed_rest;++i) {
    if (dihed_rest[i][0]==dihed[0] && dihed_rest[i][1]==dihed[1] && dihed_rest[i][2]==dihed[2] && dihed_rest[i][3]==dihed[3]) {
      flag = 1;
      break;
    }
  }

  if (flag==1) { 
    return i;
  }
  else {
    return 0;
  }
}

// 剛体中の全ての 2 面角の計算を行う関数
double pick_dihed_one_clust(int nnumaom1,
	                      int nnumaom2,
	                      int nnumaom3,
	                      int nnumaom4,
			    int nnumdihed) {
	int alpha;
	int i;
	int num_atom_n_clut;
	int num_atom_n_clut_1;
	int nnumdihedclut;
	int nnumdihedclut_1;
	int nnumdihedclut_2;
	int nnumatomorigin;
	int nnumatomtarget;

	int num_dihed_clust;

	double cn_a[3];
	double hn_a[3];
	double cn_1_a[3];
	double hn_1_a[3];

	double cc1[3];
	double cc2[3];
	double cc1cc2[3];

	double d_cn_1_a_cn_a;
	double d_hn_a_cn_a;
	double d_hn_1_a_cn_1_a ;

	double cs=0.0;
	double sn=0.0;

	double det;
//	double det[4];
	double vect[3];

	double sintheta=0.0;
	double sinphi=0.0;

	double pi;

	double theta;

//	// この結合の二面角の数を取得
//	num_dihed_clust = clust[nnumclut].f_p_clust.num_dihed_clust;

//	// 終端の剛体の場合
//	if(clust[nnumclut].terminal == terminal)
//	{
//		nnumatomorigin = clust[nnumclut].origin_atom_a-1;
//		nnumatomtarget = clust[origin_of_this_branch].terminal_atom_a[0]-1;
//	}
//	// 一つ前が終端の場合
//	else if(clust[nnumclut-1].terminal == terminal)
//	{
//		nnumatomorigin = clust[nnumclut].origin_atom_a-1;
//		nnumatomtarget = clust[origin_of_this_branch].terminal_atom_a[0]-1;
//	}
//	// 通常の剛体の場合
//	else
//	{
//		nnumatomorigin = clust[nnumclut].origin_atom_a-1;
//		nnumatomtarget = clust[nnumclut-1].terminal_atom_a[0]-1;
//	}

	// cn_a と cn_a_1 の座標の取得
	for(alpha=0;alpha<3;++alpha)
	{
		cn_a[alpha]/*a*/ = prot.coord[nnumaom2][alpha]/*a*/;
		cn_1_a[alpha]/*a*/ = prot.coord[nnumaom3][alpha]/*a*/;
	}

//	nnumdihedclut = 0;
//	for( nnumdihedclut_1 = 0;
//	     nnumdihedclut_1 < clust[nnumclut].f_p_clust.num_dihed_clust_1;
//	     ++nnumdihedclut_1 )
//	{
//		num_atom_n_clut = clust[nnumclut].f_p_clust.num_atom_n[nnumdihedclut_1]-1;

		// hn_a の座標の取得
		for(alpha=0;alpha<3;++alpha)
		{
			hn_a[alpha]/*a*/ = prot.coord[nnumaom1][alpha]/*a*/;
		}

//		for(nnumdihedclut_2 = 0; 
//		    nnumdihedclut_2 < clust[nnumclut].f_p_clust.num_dihed_clust_2; 
//		    ++nnumdihedclut_2)
//		{
			d_cn_1_a_cn_a    =0.0;
			d_hn_a_cn_a	 =0.0;
			d_hn_1_a_cn_1_a  =0.0;

			cs=0.0;
			sn=0.0;

			/*			for (alpha=0;alpha<4;++alpha)
			{
				det[alpha] = 0.0;
			}*/

		 	sintheta = 0.0;
			sinphi = 0.0;

//			num_atom_n_clut_1
//			 = clust[nnumclut].f_p_clust.num_atom_n_1[nnumdihedclut_2]-1;

			// hn_a_1 の座標の取得
			for(alpha = 0 ; alpha < 3 ; ++alpha)
			{
				hn_1_a[alpha]/*a*/
				 = prot.coord[nnumaom4][alpha]/*a*/;
			}

			for(alpha = 0 ; alpha < 3 ; ++alpha)
			{
				d_cn_1_a_cn_a/*a^2*/   +=  (cn_1_a[alpha]-cn_a[alpha])/*a*/
				                          *(cn_1_a[alpha]-cn_a[alpha])/*a*/;
				d_hn_a_cn_a/*a^2*/     +=  (hn_a[alpha]-cn_a[alpha])/*a*/
				                          *(hn_a[alpha]-cn_a[alpha])/*a*/;
				d_hn_1_a_cn_1_a/*a^2*/ +=  (hn_1_a[alpha]-cn_1_a[alpha])/*a*/
				                          *(hn_1_a[alpha]-cn_1_a[alpha])/*a*/;
			}

			d_cn_1_a_cn_a/*a*/ = sqrt(d_cn_1_a_cn_a)/*a*/;
			d_hn_a_cn_a/*a*/ = sqrt(d_hn_a_cn_a)/*a*/;
			d_hn_1_a_cn_1_a/*a*/ = sqrt(d_hn_1_a_cn_1_a)/*a*/;

			for(alpha=0;alpha<3;++alpha)
			{
				sintheta/*rad*/ += (cn_1_a[alpha]-cn_a[alpha])
				                  *(hn_a[alpha]-cn_a[alpha]);
				sinphi/*rad*/   += (cn_a[alpha]-cn_1_a[alpha])
				                  *(hn_1_a[alpha]-cn_1_a[alpha]);
			}

			sintheta/*costheta*//*rad*/ = sintheta/*rad*/
			                              /(d_cn_1_a_cn_a*d_hn_a_cn_a);
			sintheta/*rad*/ = 1.0 - sintheta/*rad*/*sintheta/*rad*/;
			if (sintheta != 0)
				sintheta/*rad*/ = sqrt(sintheta)/*rad*/;
			else 
				sintheta/*rad*/ = 0.0;

			sinphi/*rad*/ = sinphi/(d_cn_1_a_cn_a*d_hn_1_a_cn_1_a);
			sinphi/*rad*/ = 1.0 - sinphi*sinphi;
			if (sinphi != 0)
				sinphi/*rad*/ = sqrt(sinphi)/*rad*/;
			else
				sinphi/*rad*/ = 0.0;

			//			cc1[0]=(  (cn_1_a[1]-cn_a[1])
			//		 *(hn_a[2]-cn_a[2])
			//	   -  (cn_1_a[2]-cn_a[2])
			//	     *(hn_a[1]-cn_a[1]) )
			//	   /(d_cn_1_a_cn_a*d_hn_a_cn_a*sintheta);
			//
			//cc1[1]=(  (cn_1_a[2]-cn_a[2])
			//		 *(hn_a[0]-cn_a[0])
			//	   -  (cn_1_a[0]-cn_a[0])
			//	     *(hn_a[2]-cn_a[2]) )
			//	   /(d_cn_1_a_cn_a*d_hn_a_cn_a*sintheta);
			//
			//cc1[2]=(  (cn_1_a[0]-cn_a[0])
			//		 *(hn_a[1]-cn_a[1])
			//	   -  (cn_1_a[1]-cn_a[1])
			//     *(hn_a[0]-cn_a[0]) )
			//    /(d_cn_1_a_cn_a*d_hn_a_cn_a*sintheta);
			//
			//cc2[0]=(  (cn_a[1]-cn_1_a[1])
			//	     *(hn_1_a[2]-cn_1_a[2])
			//	   -  (cn_a[2]-cn_1_a[2])
			//	     *(hn_1_a[1]-cn_1_a[1]) )
			//	    /(d_cn_1_a_cn_a*d_hn_1_a_cn_1_a*sinphi);
			//
			//cc2[1]=(   (cn_a[2]-cn_1_a[2])
			//	      *(hn_1_a[0]-cn_1_a[0])
			//	   -    (cn_a[0]-cn_1_a[0])
			//	      *(hn_1_a[2]-cn_1_a[2]) )
			//	    /(d_cn_1_a_cn_a*d_hn_1_a_cn_1_a*sinphi);
			//
			//cc2[2]=(    (cn_a[0]-cn_1_a[0])
			//	       *(hn_1_a[1]-cn_1_a[1])
			//    -   (cn_a[1]-cn_1_a[1])
			//       *(hn_1_a[0]-cn_1_a[0]) )
			//    /(d_cn_1_a_cn_a*d_hn_1_a_cn_1_a*sinphi);

			cc1[0]=(  (hn_a[1]-cn_a[1])
				 *(cn_1_a[2]-cn_a[2])
				   -  
				  (hn_a[2]-cn_a[2])
				 *(cn_1_a[1]-cn_a[1])
			        )
				   /(d_cn_1_a_cn_a*d_hn_a_cn_a*sintheta);
			
			cc1[1]=(  
				   (hn_a[2]-cn_a[2])
				  *(cn_1_a[0]-cn_a[0])
				   -  
				   (hn_a[0]-cn_a[0])
				  *(cn_1_a[2]-cn_a[2])
				)
				   /(d_cn_1_a_cn_a*d_hn_a_cn_a*sintheta);
			
			cc1[2]=(  (hn_a[0]-cn_a[0])*(cn_1_a[1]-cn_a[1])
				   -(hn_a[1]-cn_a[1])*(cn_1_a[0]-cn_a[0])
				   //				    (hn_a[0]-cn_a[0])
				   //				   *(cn_1_a[1]-cn_a[1])
				   //				   -  
				   //				    (hn_a[1]-cn_a[1])
				   //				   *(cn_1_a[0]-cn_a[0])
				 )
				   /(d_cn_1_a_cn_a*d_hn_a_cn_a*sintheta);
			
			cc2[0]=(  (hn_1_a[1]-cn_1_a[1])
				    *(cn_a[2]-cn_1_a[2])
				  -  
				  (hn_1_a[2]-cn_1_a[2])
			            *(cn_a[1]-cn_1_a[1])
				)
				    /(d_cn_1_a_cn_a*d_hn_1_a_cn_1_a*sinphi);
			
			cc2[1]=(   (hn_1_a[2]-cn_1_a[2])
				      *(cn_a[0]-cn_1_a[0])
				   -    
				   (hn_1_a[0]-cn_1_a[0])
				      *(cn_a[2]-cn_1_a[2]) 
				   )
				    /(d_cn_1_a_cn_a*d_hn_1_a_cn_1_a*sinphi);
			
			cc2[2]=(    (hn_1_a[0]-cn_1_a[0])
				       *(cn_a[1]-cn_1_a[1])
				    -   
				    (hn_1_a[1]-cn_1_a[1])
				      * (cn_a[0]-cn_1_a[0]))
			    /(d_cn_1_a_cn_a*d_hn_1_a_cn_1_a*sinphi);

			for(alpha=0;alpha<3;++alpha)
			{
				cs += cc1[alpha]*cc2[alpha];
			}

			if (cs < -1.0 )
			{
				cs = -1.0;
			}
			else if (cs > 1.0 )
			{
				cs = 1.0;
			}

			cc1cc2[0] = cc1[1]*cc2[2]-cc1[2]*cc2[1];
			cc1cc2[1] = cc1[2]*cc2[0]-cc1[0]*cc2[2];
			cc1cc2[2] = cc1[0]*cc2[1]-cc1[1]*cc2[0];

			det = cc1cc2[0]*(cn_1_a[0]-cn_a[0])
				 +cc1cc2[1]*(cn_1_a[1]-cn_a[1])
				 +cc1cc2[2]*(cn_1_a[2]-cn_a[2]);

			pi = acos(-1.0);

			if (det > 0.0)
			{
			  sn=sqrt(1.0-cs*cs);
			  theta/*rad*/ = pi/*rad*/ + acos(cs)/*rad*/;
			  //	theta/*rad*/ = -acos(cs)/*rad*/;
			}
			else
			{
			  sn=-sqrt(1.0-cs*cs);
			  theta/*rad*/ = pi/*rad*/ - acos(cs)/*rad*/;
			  //	theta/*rad*/ = acos(cs)/*rad*/;
			}

//			det = cc1[0]*(hn_a[0]-hn_1_a[0])
//			     +cc1[1]*(hn_a[1]-hn_1_a[1])
// 		     	 +cc1[2]*(hn_a[2]-hn_1_a[2]);
//
//			if (det < 0)
//			{
//				theta/*rad*/ = pi/*rad*/ - acos(cs)/*rad*/;
//			}
//			else
//			{
//				theta/*rad*/ = pi/*rad*/ + acos(cs)/*rad*/;
//			}

//			det = cc1[0]*(hn_a[0]-cn_a[0])
//			     +cc1[1]*(hn_a[1]-cn_a[1])
// 		     	 +cc1[2]*(hn_a[2]-cn_a[2]);
//
//			if (det < 0)
//			{
//				theta/*rad*/ = pi/*rad*/ + acos(cs)/*rad*/;
//			}
//			else
//			{
//				theta/*rad*/ = pi/*rad*/ - acos(cs)/*rad*/;
//			}

/*			det[0]=(  (hn_1_a[1]-cn_1_a[1])
					 *(cn_a[2]-cn_1_a[2])
				   -  (hn_1_a[2]-cn_1_a[2])
				     *(cn_a[1]-cn_1_a[1]) );

			det[1]=(  (hn_1_a[2]-cn_1_a[2])
					 *(cn_a[0]-cn_1_a[0])
				   -  (hn_1_a[0]-cn_1_a[0])
				     *(cn_a[2]-cn_1_a[2]) );

			det[2]=(  (hn_1_a[0]-cn_1_a[0])
					 *(cn_a[1]-cn_1_a[1])
				   -  (hn_1_a[1]-cn_1_a[1])
				     *(cn_a[0]-cn_1_a[0]) );

			det[3] = det[0]*cn_1_a[0]
					+det[1]*cn_1_a[1]
					+det[2]*cn_1_a[2];

			for (i=0;i<3;++i)
			{
				vect[i] = cn_a[i] + cc2[i];
			}

			if (det[0]*vect[0]+det[1]*vect[1]+det[2]*vect[2]-det[3] < 0.0)
			{
				if (det[0]*hn_a[0]+det[1]*hn_a[1]+det[2]*hn_a[2]-det[3] >  0.0)
				{
					theta = pi + acos(cs);
				}
				else
				{
					theta = pi - acos(cs);
				}
			}
			else
			{
				if (det[0]*hn_a[0]+det[1]*hn_a[1]+det[2]*hn_a[2]-det[3] <  0.0)
				{
					theta = pi + acos(cs);
				}
				else
				{
					theta = pi - acos(cs);
				}
			}
*/

//	theta = acos(cs);
//			clust[nnumclut].f_p_clust.n_a_b_dihed[nnumdihedclut].dihedang/*rad*/
//			 = theta/*rad*/;

			return theta/*rad*/;


//			++nnumdihedclut;
//		}
//	}
}

// 剛体中の全ての 2 面角の計算を行う関数
double pick_dihed_one_clust4(int nnumaom1,int nnumaom2,int nnumaom3,int nnumaom4) {
  int alpha;
  
  double atom_i[3],atom_j[3],atom_k[3],atom_l[3];
  
  double vec_ij[3],vec_jk[3],vec_kl[3];
  double out_ij_jk[3],out_jk_kl[3],out_ijkl_jkkl[3];

  double d_ij_jk=0.0,d_jk_kl=0.0,d_ijkl_jkkl=0.0,d_jk=0.0;

  double det=0.0;

  double cs=0.0;
  double theta;
  double pi;

  double out_ijjk_jk[3];
  double in_ijjkjk_jkkl=0.0,in_ijjk_jkkl=0.0;

  // 座標の取得
  for(alpha=0;alpha<3;++alpha){
    atom_i[alpha] = prot.coord[nnumaom1][alpha];
    atom_j[alpha] = prot.coord[nnumaom2][alpha];
    atom_k[alpha] = prot.coord[nnumaom3][alpha];
    atom_l[alpha] = prot.coord[nnumaom4][alpha];
  }

  for (alpha=0;alpha<3;++alpha) {
    vec_ij[alpha] = atom_j[alpha]-atom_i[alpha];
    vec_jk[alpha] = atom_k[alpha]-atom_j[alpha];
    vec_kl[alpha] = atom_l[alpha]-atom_k[alpha];
  }

  out_ij_jk[0]=vec_ij[1]*vec_jk[2]-vec_ij[2]*vec_jk[1];
  out_ij_jk[1]=vec_ij[2]*vec_jk[0]-vec_ij[0]*vec_jk[2];
  out_ij_jk[2]=vec_ij[0]*vec_jk[1]-vec_ij[1]*vec_jk[0];
			
  out_jk_kl[0]=vec_jk[1]*vec_kl[2]-vec_jk[2]*vec_kl[1];
  out_jk_kl[1]=vec_jk[2]*vec_kl[0]-vec_jk[0]*vec_kl[2];
  out_jk_kl[2]=vec_jk[0]*vec_kl[1]-vec_jk[1]*vec_kl[0];

  d_ij_jk += out_ij_jk[0]*out_ij_jk[0]+out_ij_jk[1]*out_ij_jk[1]+out_ij_jk[2]*out_ij_jk[2];
  d_jk_kl += out_jk_kl[0]*out_jk_kl[0]+out_jk_kl[1]*out_jk_kl[1]+out_jk_kl[2]*out_jk_kl[2];

  d_ij_jk = sqrt(d_ij_jk);
  d_jk_kl = sqrt(d_jk_kl);

  for (alpha=0;alpha<3;++alpha) {
    out_ij_jk[alpha] = out_ij_jk[alpha]/d_ij_jk;
    out_jk_kl[alpha] = out_jk_kl[alpha]/d_jk_kl;
  }

  for(alpha=0;alpha<3;++alpha) {
    cs += out_ij_jk[alpha]*out_jk_kl[alpha];
  }

  if (cs < -1.0 ){
    cs = -1.0;
  }
  else if (cs > 1.0 ) {
    cs = 1.0;
  }

  out_ijkl_jkkl[0] = out_ij_jk[1]*out_jk_kl[2]-out_ij_jk[2]*out_jk_kl[1];
  out_ijkl_jkkl[1] = out_ij_jk[2]*out_jk_kl[0]-out_ij_jk[0]*out_jk_kl[2];
  out_ijkl_jkkl[2] = out_ij_jk[0]*out_jk_kl[1]-out_ij_jk[1]*out_jk_kl[0];

  det = out_ijkl_jkkl[0]*vec_jk[0]+out_ijkl_jkkl[1]*vec_jk[1]+out_ijkl_jkkl[2]*vec_jk[2];
  
  d_ijkl_jkkl += out_ijkl_jkkl[0]*out_ijkl_jkkl[0]+out_ijkl_jkkl[1]*out_ijkl_jkkl[1]+out_ijkl_jkkl[2]*out_ijkl_jkkl[2];
  d_ijkl_jkkl = sqrt(d_ijkl_jkkl);
  d_jk += vec_jk[0]*vec_jk[0]+vec_jk[1]*vec_jk[1]+vec_jk[2]*vec_jk[2];
  d_jk = sqrt(d_jk);

  if (d_ijkl_jkkl!=0) { /*modf 0610*/
    det = det/(d_ijkl_jkkl*d_jk);
    theta = acos(cs)*det;
  }
  else {
    theta = acos(cs);
  }
	
  pi = acos(-1.0);
  if (det<0) {
    theta = 2.0*pi+theta;
  }

  ///////////////////////////////////////////////////////////////
  /*******************************************************************/
  /* out_ijjk_jk[0] = out_ij_jk[1]*vec_jk[2]-out_ij_jk[2]*vec_jk[1]; */
  /* out_ijjk_jk[1] = out_ij_jk[2]*vec_jk[0]-out_ij_jk[0]*vec_jk[2]; */
  /* out_ijjk_jk[2] = out_ij_jk[0]*vec_jk[1]-out_ij_jk[1]*vec_jk[0]; */
  /* 								     */
  /* for (alpha=0;alpha<3;++alpha) {				     */
  /*   in_ijjkjk_jkkl += out_ijjk_jk[alpha]*out_jk_kl[alpha];	     */
  /* }								     */
  /* 								     */
  /* for (alpha=0;alpha<3;++alpha) {				     */
  /*   in_ijjk_jkkl += out_ij_jk[alpha]*out_jk_kl[alpha];	     */
  /* }								     */
  /* 								     */
  /* theta = atan2(in_ijjkjk_jkkl,in_ijjk_jkkl);		     */
  /*******************************************************************/
	
  return theta/*rad*/;

}


// 剛体中の全ての 2 面角の計算を行う関数
double pick_dihed_one_clust2(int nNumAom1,
	                      int nNumAom2,
	                      int nNumAom3,
	                      int nNumAom4,
	                      int nNumDihed)
{
	int alpha,i;

	int nnumatomorigin;
	int nnumatomtarget;
	int nnumatomtargettwo;

	double CN_A[3];
	double HN_A[3];
	double CN_1_A[3];
	double HN_1_A[3];
	double vect[3];

	double cc1[3];
	double cc2[3];

	double d1=0.0;
	double d2=0.0;
	double d4=0.0;
	double cs=0.0;

	double det[3];

	double theta;

	double sintheta=0.0;
	double sinphi=0.0;

//	// この結合の二面角の数を取得
//	num_dihed_clust = clust[nnumclut].f_p_clust.num_dihed_clust;

//	// 終端の剛体の場合
//	if(clust[nnumclut].terminal == terminal)
//	{
//		nnumatomorigin = clust[nnumclut].origin_atom_a-1;
//		nnumatomtarget = clust[origin_of_this_branch].terminal_atom_a[0]-1;
//	}
//	// 一つ前が終端の場合
//	else if(clust[nnumclut-1].terminal == terminal)
//	{
//		nnumatomorigin = clust[nnumclut].origin_atom_a-1;
//		nnumatomtarget = clust[origin_of_this_branch].terminal_atom_a[0]-1;
//	}
//	// 通常の剛体の場合
//	else
//	{
//		nnumatomorigin = clust[nnumclut].origin_atom_a-1;
//		nnumatomtarget = clust[nnumclut-1].terminal_atom_a[0]-1;
//	}

	// 原子座標の取得
	for(alpha=0;alpha<3;++alpha)
	{
		CN_A[alpha]=prot.coord[nNumAom1][alpha];
		HN_A[alpha]=prot.coord[nNumAom2][alpha];
		CN_1_A[alpha]=prot.coord[nNumAom3][alpha];
		HN_1_A[alpha]=prot.coord[nNumAom4][alpha];
	}

	// 結合長の取得_1
	for(alpha=0;alpha<3;++alpha)
	{
		d1 += (CN_1_A[alpha]-CN_A[alpha])*(CN_1_A[alpha]-CN_A[alpha]);
		d2 += (HN_A[alpha]-CN_A[alpha])*(HN_A[alpha]-CN_A[alpha]);
		d4 += (HN_1_A[alpha]-CN_1_A[alpha])*(HN_1_A[alpha]-CN_1_A[alpha]);
	}

	// 結合長の取得_2
	d1 = sqrt(d1);
	d2 = sqrt(d2);
	d4 = sqrt(d4);

	// 角(CN_1_A-CN_A-HN_A)と角(CN_1_A-HN_A-CN_A)の取得_1
	for(alpha=0;alpha<3;++alpha)
	{
		sintheta += (CN_1_A[alpha]-CN_A[alpha])*(HN_A[alpha]-CN_A[alpha]);
		sinphi += (CN_A[alpha]-CN_1_A[alpha])*(HN_1_A[alpha]-CN_1_A[alpha]);
	}

	// 角(CN_1_A-CN_A-HN_A)と角(CN_1_A-HN_A-CN_A)の取得_2
	sintheta = sintheta/(d1*d2);
	sintheta = 1.0 - sintheta*sintheta;
	if (sintheta != 0)
		sintheta = sqrt(sintheta);
	else 
		sintheta = 0.0;

	sinphi = sinphi/(d1*d4);
	sinphi = 1.0 - sinphi*sinphi;
	if (sinphi != 0)
		sinphi = sqrt(sinphi);
	else 
		sinphi = 0.0;

	// (CN_1_A-CN_A)x(HN_A-CN_A)
	cc1[0]=((CN_1_A[1]-CN_A[1])*(HN_A[2]-CN_A[2])
	       -(CN_1_A[2]-CN_A[2])*(HN_A[1]-CN_A[1]))
	       /(d1*d2*sintheta);
	cc1[1]=((CN_1_A[2]-CN_A[2])*(HN_A[0]-CN_A[0])
	       -(CN_1_A[0]-CN_A[0])*(HN_A[2]-CN_A[2]))
	       /(d1*d2*sintheta);
	cc1[2]=((CN_1_A[0]-CN_A[0])*(HN_A[1]-CN_A[1])
	       -(CN_1_A[1]-CN_A[1])*(HN_A[0]-CN_A[0]))
	       /(d1*d2*sintheta);

	// (CN_A-CN_1_A)x(HN_1_A-CN_1_A)
	cc2[0]=((CN_A[1]-CN_1_A[1])*(HN_1_A[2]-CN_1_A[2])
	       -(CN_A[2]-CN_1_A[2])*(HN_1_A[1]-CN_1_A[1]))
	       /(d1*d4*sinphi);
	cc2[1]=((CN_A[2]-CN_1_A[2])*(HN_1_A[0]-CN_1_A[0])
	       -(CN_A[0]-CN_1_A[0])*(HN_1_A[2]-CN_1_A[2]))
	       /(d1*d4*sinphi);
	cc2[2]=((CN_A[0]-CN_1_A[0])*(HN_1_A[1]-CN_1_A[1])
	       -(CN_A[1]-CN_1_A[1])*(HN_1_A[0]-CN_1_A[0]))
	       /(d1*d4*sinphi);

	for(alpha=0;alpha<3;++alpha)
	{
		cs += cc1[alpha]*cc2[alpha];
	}

	if (cs < -1.0)
	{
		cs = -1.0;
	}
	else if (cs > 1.0)
	{
		cs = 1.0;
	}

//			det = cc1[0]*(HN_A[0]-HN_1_A[0])
//			     +cc1[1]*(HN_A[1]-HN_1_A[1])
// 		     	 +cc1[2]*(HN_A[2]-HN_1_A[2]);
//
//			if (det < 0)
//			{
//				theta/*rad*/ = PI/*rad*/ - acos(cs)/*rad*/;
//			}
//			else
//			{
//				theta/*rad*/ = PI/*rad*/ + acos(cs)/*rad*/;
//			}

//			det = cc1[0]*(HN_A[0]-CN_A[0])
//			     +cc1[1]*(HN_A[1]-CN_A[1])
// 		     	 +cc1[2]*(HN_A[2]-CN_A[2]);

//			if (det < 0)
//			{
//				theta/*rad*/ = PI/*rad*/ + acos(cs)/*rad*/;
//			}
//			else
//			{
//				theta/*rad*/ = PI/*rad*/ - acos(cs)/*rad*/;
//			}


			for (i=0;i<3;++i)
			{
				vect[i] = CN_A[i] + cc2/*1*/[i];
			}

//			if (det[0]*vect[0]+det[1]*vect[1]+det[2]*vect[2]-det[3] < 0.0)
//			{
//				if (det[0]*HN_A[0]+det[1]*HN_A[1]+det[2]*HN_A[2]-det[3] > /*<*/ 0.0)
//				{
//					theta/*rad*/ = PI/*rad*/ + acos(cs)/*rad*/;
//				}
//				else
//				{
//					theta/*rad*/ = PI/*rad*/ - acos(cs)/*rad*/;
//				}
//			}
//			else
//			{
//				if (det[0]*HN_A[0]+det[1]*HN_A[1]+det[2]*HN_A[2]-det[3] < /*>*/ 0.0)
//				{
//					theta/*rad*/ = PI/*rad*/ + acos(cs)/*rad*/;
//				}
//				else
//				{
					theta/*rad*/ = PI/*rad*/ - acos(cs)/*rad*/;
//				}
//			}

	// 二面角の代入
	return theta;
/*********************************************************************************/
}

// 2 面角相互作用の計算を行う関数
void Calc_dihed_Potential_Force_for_db(void) {
  int i;
  int nNumDihed;
  int nNumClut;
  int nNumDihedType;

  double V_N_now;
  double n_dihed_now;
  double isou_dihed_now;
  double dihedang;
  double pi;

  double P_dihed_debug;
  double P_dihed_debug2=0.0;
  FILE *debug;

  potential_pro.p_dihedt = 0.0;

  for (nNumClut=0;nNumClut<prot.DOF;++nNumClut){
    potential_pro.p_dihedc[nNumClut] = 0.0;
    clust[nNumClut].f_c.f_dihed = 0.0;
  }

  for (nNumDihed=0;nNumDihed<prot.nNumDihedALL;++nNumDihed) {
    nNumDihedType = atom_dihed_pair[nNumDihed*6+4]-1;
    nNumClut = atom_dihed_pair[nNumDihed*6+5]-1;

    dihedang = pick_dihed_one_clust4(atom_dihed_pair[nNumDihed*6+0]-1,
				     atom_dihed_pair[nNumDihed*6+1]-1,
				     atom_dihed_pair[nNumDihed*6+2]-1,
				     atom_dihed_pair[nNumDihed*6+3]-1);

    dihedang = dihedang-AP.PHASE[nNumDihedType];
    pi = acos(-1.0);
    //    if (dihedang < pi) {		  
    potential_pro.p_dihedc[nNumClut] += 0.5*AP.PK[nNumDihedType]*dihedang*dihedang;
    clust[nNumClut].f_c.f_dihed += AP.PK[nNumDihedType]*4.18407*100.0*dihedang;
      //    }
    /***********************************************************************************************************/
    /* else {	  											       */
    /*   potential_pro.p_dihedc[nNumClut] += -AP.PK[nNumDihedType]*(dihedang+2.0*pi)*(dihedang+2.0*pi);	       */
    /*   clust[nNumClut].f_c.f_dihed += -2.0*AP.PK[nNumDihedType]*4.18407*100.0*(dihedang+2.0*pi);	       */
    /* }												       */
    /***********************************************************************************************************/
  }
}
