double runTACCM_CGAA_MD_NHC_MP1998_Amber_CAGo_MB(// AA /////////////////////////////////////////////////////////
						 double *crdAA,double *velAA, 
						 double *zetaAA,double *V_zetaAA, double QAA,
						 struct potential e, struct force f, double TAA, double NfKTAA,
						 double *avePEAA, double *aveKEAA,double *aveTAA,
						 double *varPEAA, double *varKEAA,double *varTAA,
						 struct my_netcdf_out_id_MCD nc_id_MCDAA,  FILE *outputfileAA,
						 // CG /////////////////////////////////////////////////////////
						 double *crdCG,double *velCG, 
						 double *zetaCG,double *V_zetaCG, double QCG,
						 struct potential_GOLM_Clementi_MB e_CG,
						 double de, double d2,
						 double TCG, double NfKTCG,
						 double *avePECG, double *aveKECG,double *aveTCG,
						 double *varPECG, double *varKECG,double *varTCG,
						 struct my_netcdf_out_id_MCD nc_id_MCDCG,  FILE *outputfileCG,
						 // Z  /////////////////////////////////////////////////////////
						 double *Z,double *velZ,double massZ,
						 double *zetaZ,double *V_zetaZ,
						 double QZ,double NfKTZ,double TZ,
						 double KZAA, double KZCG,
						 int numZ_dih,int **pairs_dihe_AA,int **pairs_dihe_CG,
						 int numZ_ang,int **pairs_angl_AA,int **pairs_angl_CG,
						 int numZ_bon,int **pairs_bond_AA,int **pairs_bond_CG,
						 double *avePEZ, double *aveKEZ,double *aveTZ,
						 double *varPEZ, double *varKEZ,double *varTZ, 
						 FILE *trjfileZ, FILE *trjfilThetaAA, FILE *trjfilThetaCG,
						 // CM  /////////////////////////////////////////////////////////
						 struct AmberParmL ap_AA,struct AmberParmL ap_CG,
						 double *mass,double *massCA,
						 int numatom, int numCAatom, 
						 int numstep, int interval,int *l,
						 double dt,double dt2,double wdt2[3] ,double wdt4[3] ,int nc,
						 double UNITT, double k_B,double pi,
						 double *PEZAA, double *PEZCG, double *PEZ) {
  int i,j,k;

  double PEAA=0.0,KEAA=0.0,EtAA,PEvAA,KEvAA;
  double PECG=0.0,KECG=0.0,EtCG,PEvCG,KEvCG;
  double KEZ,KEvZ,PEvZ,EtZ;

  double *thetaAA,*thetaCG,*frcZ,*fAA,*fCG;
  double summass,COM[3],crd_nc[MAXATOM][3];

  //  printf("nb=%5.3d na=%5.3d nd=%5.3d natom=%5.3d\n");
  //  printf("here is line 70 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");
  //  *aveKEZ=0.0;
  //  *varPEZ=0.0;
  //  *varKEZ=0.0;

  //  *aveTAA=0.0;
  //  *varTAA=0.0;

  //  *aveTCG=0.0;
  //  *varTCG=0.0;

  //  *aveTZ=0.0;
  //  *varTZ=0.0;

  ffL_calcffandforce_AA(crdAA,numatom,&e,&f);
  //  GOLM_Clementi_MB_ff_calcff(crdCG,numCAatom,de,d2,&e_CG);

  //  ffLc_calcffandforce(crdAA,numatom,&e,&f,ap_AA);
  
  fAA=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));
  fCG=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));
  frcZ=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));
  
  thetaAA=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));
  thetaCG=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));
  
  /*************************************************************/
  /* TACCM_CTheta_Amber_CAGo_MB(crdAA,numatom,thetaAA,	       */
  /* 			     numZ_dih,pairs_dihe_AA,	       */
  /* 			     numZ_ang,pairs_angl_AA,	       */
  /* 			     numZ_bon,pairs_bond_AA,	       */
  /* 			     pi);			       */
  /*************************************************************/
  /* TACCM_CTheta_Amber_CAGo_MB(crdCG,numCAatom,thetaCG,        */
  /* 			     numZ_dih,pairs_dihe_CG,	        */
  /* 			     numZ_ang,pairs_angl_CG,	        */
  /* 			     numZ_bon,pairs_bond_CG,	        */
  /* 			     pi);			        */
  /* 							        */
  /* TACCM_calc_eff_FF_Z_2(Z,numZ_dih+numZ_ang+numZ_bon,        */
  /* 			thetaAA,KZAA,fAA,pi);		        */
  /* TACCM_calc_eff_FF_Z_2(Z,numZ_dih+numZ_ang+numZ_bon,        */
  /* 			thetaCG,KZCG,fCG,pi);		        */
  /**************************************************************/

  //  TACCM_CTheta(crdAA,numatom,thetaAA,numZ_dih,pairs_dihe_AA,pi);
  /* TACCM_calc_eff_FF_Z(Z,numZ_dih,thetaAA,KZAA,fAA,pi);	    */
  /* TACCM_calc_eff_FF_Z(Z,numZ_dih,thetaCG,KZCG,fCG,pi);	    */
  /******************************************************************/
  
  //  for (i=0;i<numZ_dih+numZ_ang+numZ_bon;++i) frcZ[i]=fAA[i]+fCG[i];

  //  for (i=0;i<numZ_dih;++i) frcZ[i]=fAA[i]+fCG[i];
  
  //  printf("here is line 110 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");
  
  for (i=0;i<numstep;++i) {
    //    printf("here is line 113 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");
    //    printf("\n numstep=%5.3d\n",i);
    //    printf("nb=%5.3d na=%5.3d nd=%5.3d \n",numZ_bon,numZ_ang,numZ_dih);
    //    printf("natom=%5.3d \n",numatom);
  
    /*************************************************************/
    /* TACCM_CTheta_Amber_CAGo_MB(crdAA,numatom,thetaAA,	 */
    /* 			       numZ_dih,pairs_dihe_AA,		 */
    /* 			       numZ_ang,pairs_angl_AA,		 */
    /* 			       numZ_bon,pairs_bond_AA,		 */
    /* 			       pi);				 */
    /*************************************************************/
    //    TACCM_CTheta(crdAA,numatom,thetaAA,numZ_dih,pairs_dihe_AA,pi);
  
    //    printf("here is line 119 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");
  
    /**************************************************************/
    /* TACCM_CTheta_Amber_CAGo_MB(crdCG,numCAatom,thetaCG,	  */
    /* 			       numZ_dih,pairs_dihe_CG,		  */
    /* 			       numZ_ang,pairs_angl_CG,		  */
    /* 			       numZ_bon,pairs_bond_CG,		  */
    /* 			       pi);				  */
    /**************************************************************/
  
    //    printf("here is line 129 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");
  
    /*********************************************************************************************************/
    /* PEAA=TACCM_MD_Propagetor_NH_MP1998_AA_Amber_CAGo_MB(crdAA,velAA,mass,zetaAA,V_zetaAA,		     */
    /* 							QAA,NfKTAA,numatom,&KEAA,&KEvAA,&PEvAA,		     */
    /* 							dt,dt2,nc,wdt4,wdt2,				     */
    /* 							&e,&f,Z,numZ_dih,numZ_ang,numZ_bon,		     */
    /* 							thetaAA,KZAA,					     */
    /* 							pairs_dihe_AA,pairs_angl_AA,pairs_bond_AA,	     */
    /* 							PEZAA,pi);					     */
    /*********************************************************************************************************/

   /***********************************************************************************************/
   /* PEAA=TACCM_MDc_Propagetor_NH_MP1998_AA_test(crdAA,velAA,mass,zetaAA,V_zetaAA,		  */
   /* 						QAA,NfKTAA,numatom,&KEAA,&KEvAA,&PEvAA,		  */
   /* 						dt,dt2,nc,wdt4,wdt2,				  */
   /* 						&e,&f,						  */
   /* 						ap_AA,Z,numZ_dih,				  */
   /* 						thetaAA,KZAA,pairs_dihe_AA,PEZAA,pi);		  */
   /***********************************************************************************************/
  
    //    printf("here is line 138 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");
  
    /*********************************************************************************************************/
    /* PECG=TACCM_MD_Propagetor_NH_MP1998_CG_Amber_CAGo_MB(crdCG,velCG,massCA,zetaCG,V_zetaCG,		     */
    /* 							QCG,NfKTCG,numCAatom,&KECG,&KEvCG,&PEvCG,	     */
    /* 							dt,dt2,nc,wdt4,wdt2,				     */
    /* 							&e_CG,de,d2,Z,numZ_dih,numZ_ang,numZ_bon,	     */
    /* 							thetaCG,KZCG,					     */
    /* 							pairs_dihe_CG,pairs_angl_CG,pairs_bond_CG,	     */
    /* 							PEZCG,pi);					     */
    /*********************************************************************************************************/
  
    //    printf("here is line 147 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");
  
    /*********************************************************************************************/
    /* TACCM_CGAA_MD_Propagetor_NH_MP1998_Z_2(Z,velZ,massZ,thetaAA,thetaCG,zetaZ,V_zetaZ,	 */
    /* 					   QZ,NfKTZ,numZ_dih+numZ_ang+numZ_bon,			 */
    /* 					   &KEZ,&KEvZ,&PEvZ,					 */
    /* 					   dt,dt2,nc,wdt4,wdt2,KZAA,KZCG,			 */
    /* 					   PEZAA,PEZCG,PEZ,frcZ,pi);				 */
    /*********************************************************************************************/
  
    //    printf("here is line 156 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");
  
    if (i%interval==0) {
  
      //      printf("here is line 165 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");
  
      KEAA=KEAA/UNITT;     TAA=KEAA/((3*numatom)*k_B)*2.0;
      PEvAA=PEvAA/UNITT;   KEvAA=KEvAA/UNITT;
  
      KECG=KECG/UNITT;     TCG=KECG/((3*numCAatom)*k_B)*2.0;
      PEvCG=PEvCG/UNITT;  KEvCG=KEvCG/UNITT;
  
      PEAA=e.p_d_t+e.p_a_t+e.p_b_t;
      EtAA=PEAA+KEAA+PEvAA+KEvAA;
      fprintf(outputfileAA,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "
    	                    ,i+1,PEAA,KEAA,KEvAA,PEvAA,EtAA,TAA);
  
      PECG=e_CG.p_MB;
      EtCG=PECG+KECG+PEvCG+KEvCG;
      fprintf(outputfileCG,"%d %8.4e %8.4e %8.4e %8.4e %8.4e %8.4e "
    	                  ,i+1,PECG,KECG,KEvCG,PEvCG,EtCG,TCG);
  
      //      *avePEAA=(i*(*avePEAA)+PEAA)/(i+1); *varPEAA=(i*(*varPEAA)+PEAA*PEAA)/(i+1);
      //      *aveKEAA=(i*(*aveKEAA)+KEAA)/(i+1); *varKEAA=(i*(*varKEAA)+KEAA*KEAA)/(i+1);
      //      *aveTAA=(i*(*aveTAA)+TAA)/(i+1);  *varTAA=(i*(*varTAA)+TAA*TAA)/(i+1);
  
      //      *avePECG=(i*(*avePECG)+PECG)/(i+1); *varPECG=(i*(*varPECG)+PECG*PECG)/(i+1);
      //      *aveKECG=(i*(*aveKECG)+KECG)/(i+1); *varKECG=(i*(*varKECG)+KECG*KECG)/(i+1);
      //      *aveTCG=(i*(*aveTCG)+TCG)/(i+1);  *varTCG=(i*(*varTCG)+TCG*TCG)/(i+1);
  
      summass=0.0; for (j=0;j<numatom;++j) summass+=mass[j];
      for (j=0;j<3;++j) COM[j]=0.0;
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) COM[k]+=mass[j]*crdAA[j*3+k]/summass;
      for (j=0;j<numatom;++j)  for (k=0;k<3;++k) crdAA[j*3+k]-=COM[k];
      for (j=0;j<numatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdAA[j*3+k];
      myncL_put_crd_MCD(nc_id_MCDAA,*l,crd_nc);
  
      summass=0.0; for (j=0;j<numCAatom;++j) summass+=massCA[j];
      for (j=0;j<3;++j) COM[j]=0.0;
      for (j=0;j<numCAatom;++j)  for (k=0;k<3;++k) COM[k]+=massCA[j]*crdCG[j*3+k]/summass;
      for (j=0;j<numCAatom;++j)  for (k=0;k<3;++k) crdCG[j*3+k]-=COM[k];
      for (j=0;j<numCAatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdCG[j*3+k];
      for (j=0;j<numCAatom;++j) for (k=0;k<3;++k) crd_nc[j][k]=crdCG[j*3+k];
      myncL_put_crd_MCD(nc_id_MCDCG,*l,crd_nc);
      ++(*l);
  
      ///////////////// TACCM //////////////////////
      KEZ=KEZ/UNITT;      TZ=KEZ/((numZ_dih+numZ_ang+numZ_bon)*k_B)*2.0;
  
      PEvZ=PEvZ/UNITT;      KEvZ=KEvZ/UNITT;
  
      EtZ=*PEZ+KEZ+PEvZ+KEvZ;
      fprintf(outputfileAA,"%d %e %e %e %e %e %e %e\n",i+1,*PEZ,KEZ,KEvZ,PEvZ,EtZ,TZ);
  
      for (j=0;j<(numZ_dih+numZ_ang+numZ_bon);++j) fprintf(trjfileZ,"%e ",Z[j]);
      fprintf(trjfileZ,"\n");
      for (j=0;j<(numZ_dih+numZ_ang+numZ_bon);++j) fprintf(trjfilThetaAA,"%e ",thetaAA[j]);
      fprintf(trjfilThetaAA,"\n");
      for (j=0;j<(numZ_dih+numZ_ang+numZ_bon);++j) fprintf(trjfilThetaCG,"%e ",thetaCG[j]);
      fprintf(trjfilThetaCG,"\n");
      ///////////////// TACCM //////////////////////
    }
  }
  //  printf("here is line 246 on TACCM_CGAA_MDrun_Amber_CAGo_MB.c\n");

  return *PEZ;
}

double TACCM_MD_Propagetor_NH_MP1998_AA_Amber_CAGo_MB(double *crd,double *vel,double *mass,
						      double *zeta,double *V_zeta,double Q,
						      double NfKT,int numatom,double *KE,double *KEv,double *PEv,
						      double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
						      struct potential *e, struct force *f,
						      double *Z, int numZ_dih,int numZ_ang, int numZ_bon,
						      double *theta,double Kapa,
						      int **pairs_dihe_AA, int **pairs_angl_AA, int **pairs_bond_AA,
						      double *PEZ,double pi) {
  int i,j,k;
  double *frc;
  double **frcZ;

  frc=(double *)gcemalloc(sizeof(double)*numatom*3);
  frcZ=(double **)gcemalloc(sizeof(double *)*numatom);
  for (i=0;i<numatom;++i) frcZ[i]=(double *)gcemalloc(sizeof(double)*3);

  /*************************************************************************************************/
  /* TACCM_CTheta_Amber_CAGo_MB(crd,numatom,theta,						   */
  /* 		       numZ_dih,pairs_dihe_AA,							   */
  /* 		       numZ_ang,pairs_angl_AA,							   */
  /* 		       numZ_bon,pairs_bond_AA,							   */
  /* 		       pi);									   */
  /* *PEZ=TACCM_calc_eff_FF_MD_Amber_CAGo_MB(crd,numatom,theta,Z,				   */
  /* 					  numZ_dih,numZ_ang,numZ_bon,				   */
  /* 					  Kapa,frcZ,						   */
  /* 					  pairs_dihe_AA,pairs_angl_AA,pairs_bond_AA,pi);	   */
  /*************************************************************************************************/

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]/*+frcZ[i][j]*/;

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];

  ffL_calcffandforce_AA(crd,numatom,e,f);
  /*************************************************************************************************/
  /* TACCM_CTheta_Amber_CAGo_MB(crd,numatom,theta,						   */
  /* 		       numZ_dih,pairs_dihe_AA,							   */
  /* 		       numZ_ang,pairs_angl_AA,							   */
  /* 		       numZ_bon,pairs_bond_AA,							   */
  /* 		       pi);									   */
  /* *PEZ=TACCM_calc_eff_FF_MD_Amber_CAGo_MB(crd,numatom,theta,Z,				   */
  /* 					  numZ_dih,numZ_ang,numZ_bon,				   */
  /* 					  Kapa,frcZ,						   */
  /* 					  pairs_dihe_AA,pairs_angl_AA,pairs_bond_AA,pi);	   */
  /*************************************************************************************************/

  for (i=0;i<numatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=-(*f).f_b[i*3+j]+(*f).f_a[i*3+j]+(*f).f_d[i*3+j]/*+frcZ[i][j]*/;
  for (i=0;i<numatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numatom,nc,wdt4,wdt2);

  *KE=0.0; for (i=0;i<numatom;++i) for (j=0;j<3;++j) *KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return 0.0;
}

double TACCM_MD_Propagetor_NH_MP1998_CG_Amber_CAGo_MB(double *crd,double *vel,double *mass,
						      double *zeta,double *V_zeta,double Q,
						      double NfKT,int numCAatom,double *KE,double *KEv,double *PEv,
						      double dt,double dt2,int nc,double wdt4[3],double wdt2[3],
						      struct potential_GOLM_Clementi_MB *e_CG,
						      double de, double d2,
						      double *Z, int numZ_dih,int numZ_ang, int numZ_bon,
						      double *theta,double Kapa,
						      int **pairs_dihe_CG, int **pairs_angl_CG, int **pairs_bond_CG,
						      double *PEZ,double pi) {
  int i,j,k;
  double *frc;
  double **frcZ;

  frc=(double *)gcemalloc(sizeof(double)*numCAatom*3);
  frcZ=(double **)gcemalloc(sizeof(double *)*numCAatom);
  for (i=0;i<numCAatom;++i) frcZ[i]=(double *)gcemalloc(sizeof(double)*3);

  TACCM_CTheta_Amber_CAGo_MB(crd,numCAatom,theta,
		       numZ_dih,pairs_dihe_CG,
		       numZ_ang,pairs_angl_CG,
		       numZ_bon,pairs_bond_CG,
		       pi);
  *PEZ=TACCM_calc_eff_FF_MD_Amber_CAGo_MB(crd,numCAatom,theta,Z,
					  numZ_dih,numZ_ang,numZ_bon,
					  Kapa,frcZ,
					  pairs_dihe_CG,pairs_angl_CG,pairs_bond_CG,pi);

  for (i=0;i<numCAatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_CG).f_MB[i][j]+frcZ[i][j];

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numCAatom,nc,wdt4,wdt2);

  //////////////////////////////////////////////////////////////////////////////////////
  for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) crd[i*3+j]+=dt*vel[i*3+j];

  GOLM_Clementi_MB_ff_calcff(crd,numCAatom,de,d2,e_CG);
  TACCM_CTheta_Amber_CAGo_MB(crd,numCAatom,theta,
		       numZ_dih,pairs_dihe_CG,
		       numZ_ang,pairs_angl_CG,
		       numZ_bon,pairs_bond_CG,
		       pi);
  *PEZ=TACCM_calc_eff_FF_MD_Amber_CAGo_MB(crd,numCAatom,theta,Z,
					  numZ_dih,numZ_ang,numZ_bon,
					  Kapa,frcZ,
					  pairs_dihe_CG,pairs_angl_CG,pairs_bond_CG,pi);

  for (i=0;i<numCAatom;++i)
    for (j=0;j<3;++j) 
      frc[i*3+j]=(*e_CG).f_MB[i][j]+frcZ[i][j];
  for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) vel[i*3+j]+=dt2*frc[i*3+j]/mass[i];
  //////////////////////////////////////////////////////////////////////////////////////

  MD_Propagetor_NH_Single_part_MP1996(vel,mass,zeta,V_zeta,Q,NfKT,numCAatom,nc,wdt4,wdt2);

  *KE=0.0; for (i=0;i<numCAatom;++i) for (j=0;j<3;++j) *KE+=0.5*mass[i]*vel[i*3+j]*vel[i*3+j];

  *KEv=0.5*Q*(*V_zeta)*(*V_zeta);
  *PEv=NfKT*(*zeta);

  return 0.0;
}

double CE_TACCMb_CGAA_Amber_CAGo_MB(double *crdAA,double *crdCG,double *Z, int numatom,int numCAatom,
				    //int numZ_dih, int **pairs_dihed_AA, int **pairs_dihed_CG,
				    //int numZ_ang, int **pairs_angle_AA, int **pairs_angle_CG,
				    //int numZ_bon, int **pairs_bond_AA,  int **pairs_bond_CG,
				    int numZ_dih, int numZ_ang, int numZ_bon, 
				    int **pairs_dihed_AA, int **pairs_dihed_CG,
				    int **pairs_angle_AA, int **pairs_angle_CG,
				    int **pairs_bond_AA,  int **pairs_bond_CG,
				    double KZAA,double KZCG, double pi,
				    double *EAA,double *ECG,double *EZ) {
  int i,j;
  //  int numZ;
  double *thetaAA,*thetaCG;
  double delta;

  printf("numZ_bon=%2d\n",numZ_bon);
  printf("numZ_ang=%2d\n",numZ_ang);
  printf("numZ_dih=%2d\n",numZ_dih);

  printf("test1:%4d -%4d  @TACCM_CGAAMDrun_Amber_CAGo_MB.c\n"
  	 ,pairs_bond_AA[0][0],pairs_bond_AA[0][1]);

  printf("test2:%4d -%4d  @TACCM_CGAAMDrun_Amber_CAGo_MB.c\n"
  	 ,pairs_bond_AA[1][0],pairs_bond_AA[1][1]);

  printf("test3:%4d -%4d  @TACCM_CGAAMDrun_Amber_CAGo_MB.c\n"
  	 ,pairs_bond_AA[2][0],pairs_bond_AA[2][1]);

  printf("test1:%4d -%4d  @TACCM_CGAAMDrun_Amber_CAGo_MB.c\n"
  	 ,pairs_bond_CG[0][0],pairs_bond_CG[0][1]);

  printf("test2:%4d -%4d  @TACCM_CGAAMDrun_Amber_CAGo_MB.c\n"
  	 ,pairs_bond_CG[1][0],pairs_bond_CG[1][1]);

  printf("test3:%4d -%4d  @TACCM_CGAAMDrun_Amber_CAGo_MB.c\n"
  	 ,pairs_bond_CG[2][0],pairs_bond_CG[2][1]);

  thetaAA=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));
  thetaCG=(double *)gcemalloc(sizeof(double)*(numZ_dih+numZ_ang+numZ_bon));
  //  printf("here is line 614 on TACCM_CGAAMDrun_Amber_CAGo_MB.c\n");
  //  printf("%5.3d %5.3d %5.3d\n",numZ_bon,numZ_ang,numZ_dih);
  /*****************************************************************/
  /* TACCM_CTheta_Amber_CAGo_MB(crdAA,numatom,thetaAA,		   */
  /* 			     numZ_dih,pairs_dihed_AA,		   */
  /* 			     numZ_ang,pairs_angle_AA,		   */
  /* 			     numZ_bon,pairs_bond_AA,pi);	   */
  /*****************************************************************/

  printf("test1:%4d -%4d -%4d @TACCM_CGAAMDrun_Amber_CAGo_MB.c\n"
  	 ,pairs_angle_CG[0][0],pairs_angle_CG[0][1],pairs_angle_CG[0][2]);

  printf("test2:%4d -%4d -%4d @TACCM_CGAAMDrun_Amber_CAGo_MB.c\n"
  	 ,pairs_angle_CG[1][0],pairs_angle_CG[1][1],pairs_angle_CG[1][2]);

  printf("test3:%4d -%4d -%4d @TACCM_CGAAMDrun_Amber_CAGo_MB.c\n"
  	 ,pairs_angle_CG[2][0],pairs_angle_CG[2][1],pairs_angle_CG[2][2]);

  printf("test4:%4d -%4d -%4d @TACCM_CGAAMDrun_Amber_CAGo_MB.c\n"
  	 ,pairs_angle_CG[3][0],pairs_angle_CG[3][1],pairs_angle_CG[3][2]);

  /***************************************************************************************************************/
  /* for (i=0;i<numZ_dih;++i) {											 */
  /*   printf("%4d -%4d -%4d -%4d @TACCM_CGAAMDrun_Amber_CAGo_MB.c\n"						 */
  /* 	   ,pairs_dihed_CG[i][0],pairs_dihed_CG[i][1]								 */
  /* 	   ,pairs_dihed_CG[i][2],pairs_dihed_CG[i][3]);								 */
  /* }														 */
  /* 														 */
  /* printf("\n");												 */
  /* 														 */
  /* for (i=0;i</\*1*\/numZ_ang;++i) {										 */
  /*   printf("A%d:%4d -%4d -%4d @TACCM_CGAAMDrun_Amber_CAGo_MB.c\n"						 */
  /* 	   ,i,pairs_angle_CG[i][0],pairs_angle_CG[i][1],pairs_angle_CG[i][2]);					 */
  /* }														 */
  /* 														 */
  /* printf("\n");												 */
  /* 														 */
  /* for (i=0;i</\*1*\/numZ_bon;++i) {										 */
  /*   printf("B%d:%4d -%4d @TACCM_CGAAMDrun_Amber_CAGo_MB.c\n",pairs_bond_CG[i][0],pairs_bond_CG[i][1]);	 */
  /* }														 */
  /***************************************************************************************************************/
    
  TACCM_CTheta_Amber_CAGo_MB(crdCG,numCAatom,thetaCG,
  			     numZ_dih,pairs_dihed_CG,
  			     numZ_ang,pairs_angle_CG,
  			     numZ_bon,pairs_bond_CG,pi);
  //  printf("here is line 623 on TACCM_CGAAMDrun_Amber_CAGo_MB.c\n");
  *EAA=0.0;
  for (i=0;i<numZ_bon;++i) {
    delta=Z[i]-thetaAA[i];
    *EAA+=0.5*KZAA*delta*delta;
  }
  for (i=0;i<numZ_ang;++i) {
    if ((delta=Z[i+numZ_bon]-thetaAA[i+numZ_bon])>pi) delta-=2.0*pi;
    else if ((delta=Z[i+numZ_bon]-thetaAA[i+numZ_bon])<-1.0*pi) delta+=2.0*pi;
    *EAA+=0.5*KZAA*delta*delta;
  }
  for (i=0;i<numZ_dih;++i) {
    if ((delta=Z[i+numZ_bon+numZ_ang]-thetaAA[i+numZ_bon+numZ_ang])>pi) delta-=2.0*pi;
    else if ((delta=Z[i+numZ_bon+numZ_ang]-thetaAA[i+numZ_bon+numZ_ang])<-1.0*pi) delta+=2.0*pi;
    *EAA+=0.5*KZAA*delta*delta;
  }
  
  *ECG=0.0;
  for (i=0;i<numZ_bon;++i) {
    delta=Z[i]-thetaCG[i];
    *ECG+=0.5*KZCG*delta*delta;
  }
  for (i=0;i<numZ_ang;++i) {
    if ((delta=Z[i+numZ_bon]-thetaCG[i+numZ_bon])>pi) delta-=2.0*pi;
    else if ((delta=Z[i+numZ_bon]-thetaCG[i+numZ_bon])<-1.0*pi) delta+=2.0*pi;
    *ECG+=0.5*KZCG*delta*delta;
  }
  for (i=0;i<numZ_dih;++i) {
    if ((delta=Z[i+numZ_bon+numZ_ang]-thetaCG[i+numZ_bon+numZ_ang])>pi) delta-=2.0*pi;
    else if ((delta=Z[i+numZ_bon+numZ_ang]-thetaCG[i+numZ_bon+numZ_ang])<-1.0*pi) delta+=2.0*pi;
    *ECG+=0.5*KZCG*delta*delta;
  }
  
  *EZ=(*EAA)+(*ECG);
  //  printf("here is line 657 on TACCM_CGAAMDrun_Amber_CAGo_MB.c\n");
}
