
double **MPI_CGAATREM_TACCM_MD_pep_NHC_MP1998_Amber_CAGo_MB(int myrank,int num_procs,int tag, MPI_Status* status,
							    int numRE, int numEX, double *KZAA, double *KZCG,
							    struct AADataforREMD_Amber_CAGo_MB AAData,
							    struct CGDataforREMD_Amber_CAGo_MB CGData,
							    struct TACCMDataforREMD_Amber_CAGo_MB ZData,
							    struct AACGCommonDataforREMD_Amber_CAGo_MB CData,
							    struct AmberParmL ap_AA, struct AmberParmL ap_CG,
							    double de, double d2,
							    double T0AA,double T0CG, double T0Z, 
							    int numstep, int interval, 
							    double dt,double dt2,
							    double wdt2[3],double wdt4[3], int nc,
							    double UNITT, double k_B, double tau, double pi,
							    double parameterCG, FILE* logfile) {
  int i,j,k,c,l=0;
  int m,n;
  int your_rank;
  int count;

  double betaAA,betaCG,betaZ;

  double EAAm_Xi,EAAm_Xj,EAAn_Xi,EAAn_Xj;
  double ECGm_Xi,ECGm_Xj,ECGn_Xi,ECGn_Xj;
  double EZm_Xi,EZm_Xj,EZn_Xi,EZn_Xj;

  double KZAAn,KZCGn;

  double *crdAA_receiv,*crdCG_receiv;
  double *Z_receiv;
  double *crd_temp,*vel_temp;
  double *Z_temp,*vel_Z_temp;

  double *SData,*RData,*SKZ,*RKZ;
  double delta=0.0;

  int exchange_data[4],exchange_data_receiv[4],*exchange_data_sum,*exchanged_data,*exchanged_data_receiv;
  double **acc_ratio;

  unsigned long init[4]={0x123, 0x234, 0x345, 0x456},length=4;
  init_by_array(init,length);

  betaAA=1.0/(k_B*T0AA);
  betaCG=1.0/(k_B*T0CG);
  betaZ=1.0/(k_B*T0Z);

  SData=(double *)gcemalloc(sizeof(double)*6/*8*/);
  RData=(double *)gcemalloc(sizeof(double)*6/*8*/);
  SKZ=(double *)gcemalloc(sizeof(double)*2);
  RKZ=(double *)gcemalloc(sizeof(double)*2);

  exchange_data_sum=(int *)gcemalloc(sizeof(int)*numRE*4);
  exchanged_data=(int *)gcemalloc(sizeof(int)*numRE*2);
  exchanged_data_receiv=(int *)gcemalloc(sizeof(int)*numRE*2);

  acc_ratio=(double **)gcemalloc(sizeof(double *)*numRE);
  for (i=0;i<numRE;++i) acc_ratio[i]=(double *)gcemalloc(sizeof(double)*numRE);
  for (i=0;i<numRE;++i) for (j=0;j<numRE;++j) acc_ratio[i][j]=0.0;

  //  printf("here is line 74 on REMDCGAA_TACCM_MPI_2_Amber-CAGo-MB.c\n");

  //  printf("nd=%5.3d na=%5.3d nb=%5.3d\n",ZData.numZ_dih,ZData.numZ_ang,ZData.numZ_bon);


  for (i=0;i<numEX;++i) {
        runTACCM_CGAA_MD_NHC_MP1998_Amber_CAGo_MB(// AA /////////////////////////////////////////////////////////
						  AAData.crd,AAData.vel,
						  &(AAData.zeta),&(AAData.V_zeta),AAData.Q,
						  AAData.e,AAData.f,AAData.T,AAData.NfKT,
						  AAData.avePE,AAData.aveKE,AAData.aveT,
						  AAData.varPE,AAData.varKE,AAData.varT,
						  AAData.nc_id_MCD,AAData.outputfile,
						  // CG /////////////////////////////////////////////////////////
						  CGData.crd,CGData.vel,
						  &(CGData.zeta),&(CGData.V_zeta),CGData.Q,
						  CGData.e_GOLM,
						  de,d2,
						  CGData.T,CGData.NfKT,
						  CGData.avePE,CGData.aveKE,CGData.aveT,
						  CGData.varPE,CGData.varKE,CGData.varT,
						  CGData.nc_id_MCD,CGData.outputfile,
						  // Z  /////////////////////////////////////////////////////////
						  ZData.Z,ZData.velZ,ZData.massZ,
						  &(ZData.zetaZ),&(ZData.V_zetaZ),
						  ZData.QZ,ZData.NfKTZ,ZData.T,
						  ZData.KZAA,ZData.KZCG,
						  ZData.numZ_dih,ZData.pairs_dihed_AA,ZData.pairs_dihed_CG,
						  ZData.numZ_ang,ZData.pairs_angle_AA,ZData.pairs_angle_CG,
						  ZData.numZ_bon,ZData.pairs_bond_AA,ZData.pairs_bond_CG,
						  ZData.avePEZ,ZData.aveKEZ,ZData.aveTZ,
						  ZData.varPEZ,ZData.varKEZ,ZData.varTZ,
						  ZData.trjfileZ,ZData.trjfilThetaAA,ZData.trjfilThetaCG,
						  // CM  /////////////////////////////////////////////////////////
						  ap_AA,ap_CG,
						  CData.mass,CData.massCA,
						  (CData.numatom),(CData.numCAatom),
						  numstep,interval,&l,
						  dt,dt2,wdt2,wdt4,nc,UNITT,k_B,pi,
						  &EAAm_Xi,&ECGm_Xi,&EZm_Xi);

    //    printf("here is line 112 on REMDCGAA_TACCM_MPI_2_Amber-CAGo-MB.c\n");
    m=REMD_purmutation_func(myrank);
    //    printf("here is line 116 on REMDCGAA_TACCM_MPI_2_Amber-CAGo-MB.c\n");
    fprintf(logfile,"%d-th: m= %d KZAA= %8.4lf KZCG= %8.4lf\n",i+1,m,ZData.KZAA,ZData.KZCG);
    //    printf("here is line 118 on REMDCGAA_TACCM_MPI_2_Amber-CAGo-MB.c\n");
    
    AAData.zeta=0.0; AAData.V_zeta=0.0;
    CGData.zeta=0.0; CGData.V_zeta=0.0;
    ZData.zetaZ=0.0; ZData.V_zetaZ=0.0;
    
    if(m%2==0 ) {
      if ( i%2 == 0 ) n=m+1;
      else {
    	if ( m==0 )  n=numRE-1;
    	else  n=m-1;
      }
    }
    else {
      if ( i%2 == 0 ) n=m-1;
      else {
    	if ( m==numRE-1 )  n=0;
    	else  n=m+1;
      }
    }
    
    //    printf("here is line 136 on REMDCGAA_TACCM_MPI_2_Amber-CAGo-MB.c\n");
    your_rank=REMD_purmutation_inverse_func(n);
    //    printf("here is line 138 on REMDCGAA_TACCM_MPI_2_Amber-CAGo-MB.c\n");
    
    SKZ[0]=ZData.KZAA;
    SKZ[1]=ZData.KZCG;
    //    printf("107\n");
    MPI_Send(SKZ, 2, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    MPI_Recv(RKZ, 2, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);
    //    printf("110\n");
    KZAAn=RKZ[0];
    KZCGn=RKZ[1];
    //    printf("here is line 147 on REMDCGAA_TACCM_MPI_2_Amber-CAGo-MB.c\n");
    MPI_Barrier(MPI_COMM_WORLD);
    //    printf("here is line 149 on REMDCGAA_TACCM_MPI_2_Amber-CAGo-MB.c\n");
    //    printf("%5.3d %5.3d %5.3d\n",ZData.numZ_bon,ZData.numZ_ang,ZData.numZ_dih);

    /*********************************************************************************************************/
    /* for (i=0;i<ZData.numZ_dih;++i) {									     */
    /*   printf("D%d:%4d -%4d -%4d -%4d REMDCGAA_TACCM_MPI_2_Amber-CAGo-MB.c\n"				     */
    /* 	     ,i,ZData.pairs_dihed_CG[i][0],ZData.pairs_dihed_CG[i][1]					     */
    /* 	     ,ZData.pairs_dihed_CG[i][2],ZData.pairs_dihed_CG[i][3]);					     */
    /* }												     */
    /* 													     */
    /* printf("\n");											     */
    /* for (i=0;i<ZData.numZ_bon;++i) {									     */
    /*   printf("B%d:%4d -%4d @REMDCGAA_TACCM_MPI_2_Amber-CAGo-MB.c\n"					     */
    /* 	     ,i,ZData.pairs_bond_CG[i][0],ZData.pairs_bond_CG[i][1]);					     */
    /* }												     */
    /* 													     */
    /* printf("\n");											     */
    /* for (i=0;i<ZData.numZ_ang;++i) {									     */
    /*   printf("A%d:%4d -%4d -%4d REMDCGAA_TACCM_MPI_2_Amber-CAGo-MB.c\n"				     */
    /* 	     ,i,ZData.pairs_angle_CG[i][0],ZData.pairs_angle_CG[i][1],ZData.pairs_angle_CG[i][2]);	     */
    /* }												     */
    /*********************************************************************************************************/
    
    /********************************************************************************************************/
    /* CE_TACCMb_CGAA_Amber_CAGo_MB(AAData.crd,CGData.crd,ZData.Z,(CData.numatom),(CData.numCAatom),	    */
    /* 				 //ZData.numZ_dih,ZData.pairs_dihed_AA,ZData.pairs_dihed_CG,		    */
    /* 				 //ZData.numZ_ang,ZData.pairs_dihed_AA,ZData.pairs_dihed_CG,		    */
    /* 				 //ZData.numZ_bon,ZData.pairs_dihed_AA,ZData.pairs_dihed_CG,		    */
    /* 				 ZData.numZ_dih,ZData.numZ_ang,ZData.numZ_bon,				    */
    /* 				 ZData.pairs_dihed_AA,ZData.pairs_dihed_CG,				    */
    /* 				 ZData.pairs_dihed_AA,ZData.pairs_dihed_CG,				    */
    /* 				 ZData.pairs_dihed_AA,ZData.pairs_dihed_CG,				    */
    /* 				 KZAAn,KZCGn,pi,							    */
    /* 				 &EAAn_Xi,&ECGn_Xi,&EZn_Xi);						    */
    /********************************************************************************************************/
    //    printf("here is line 158 on REMDCGAA_TACCM_MPI_2_Amber-CAGo-MB.c\n");
    SData[0]=EAAm_Xi;
    SData[1]=ECGm_Xi;
    SData[2]=EZm_Xi;
    
    SData[3]=EAAn_Xi;
    SData[4]=ECGn_Xi;
    SData[5]=EZn_Xi;
    
    //    SData[6]=ZData.KZAA;
    //    SData[7]=ZData.KZCG;
    //    printf("148\n");
    MPI_Send(SData, /*8*/6, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    MPI_Recv(RData, /*8*/6, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);
    //    printf("151\n");
    EAAn_Xj=RData[0];
    ECGn_Xj=RData[1];
    EZn_Xj=RData[2];
    
    EAAm_Xj=RData[3];
    ECGm_Xj=RData[4];
    EZm_Xj=RData[5];
    
    //    KZAAn=RData[6];
    //    KZCGn=RData[7];
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (n>m) {
      //      delta=betaAA*((EAAm_Xj-EAAm_Xi)-(EAAn_Xi-EAAn_Xj))
      //      	+betaCG*((ECGm_Xj-ECGm_Xi)-(ECGn_Xi-ECGn_Xj))
      //      	+betaZ*((EZm_Xj-EZm_Xi)-(EZn_Xi-EZn_Xj));
    
      delta=betaAA*((EAAm_Xj-EAAm_Xi)-(EAAn_Xj-EAAn_Xi))
      	+betaCG*((ECGm_Xj-ECGm_Xi)-(ECGn_Xj-ECGn_Xi))
      	+betaZ*((EZm_Xj-EZm_Xi)-(EZn_Xj-EZn_Xi));
    
      /**************************************************************/
      /* delta=-betaAA*((EAAm_Xj-EAAm_Xi)-(EAAn_Xi-EAAn_Xj))	    */
      /* 	-betaCG*((ECGm_Xj-ECGm_Xi)-(ECGn_Xi-ECGn_Xj))	    */
      /* 	-betaZ*((EZm_Xj-EZm_Xi)-(EZn_Xi-EZn_Xj));	    */
      /**************************************************************/
    
      MPI_Send(&delta, 1, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD);
    }
    else MPI_Recv(&delta, 1, MPI_DOUBLE, your_rank, tag, MPI_COMM_WORLD, status);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if((c=Metropolis(delta))==1) {
      if (myrank<your_rank)
    	printf("%d-th Acc %d-%d : /_\\ = %8.4lf \n EAAmi=%8.4lf ECGmi=%8.4lf EZmi=%8.4lf EAAni=%8.4lf ECGni=%8.4lf EZni=%8.4lf\n EAAmj=%8.4lf ECGmj=%8.4lf EZmj=%8.4lf EAAnj=%8.4lf ECGnj=%8.4lf EZnj=%8.4lf\n KZAAm=%8.4lf KZAAn=%8.4lf KZCGm=%8.4lf KZCGn=%8.4lf\n",
    	       i+1,myrank,your_rank,myrank,delta,
    	       EAAm_Xi,ECGm_Xi,EZm_Xi,EAAn_Xi,ECGn_Xi,EZn_Xi,
    	       EAAm_Xj,ECGm_Xj,EZm_Xj,EAAn_Xj,ECGn_Xj,EZn_Xj,
    	       ZData.KZAA,ZData.KZCG,KZAAn,KZCGn);
      ZData.KZAA=KZAAn;
      ZData.KZCG=KZCGn;
    
      if (myrank!=0) {
    	exchange_data[0]=myrank;
    	exchange_data[1]=your_rank;
    	exchange_data[2]=m;
    	exchange_data[3]=n;
      }
      else {
    	REMD_exchange_purmutation_funcs(myrank,your_rank,m,n);
    	//	printf("176 m=%d n=%d\n",m,n);
    	acc_ratio[m][n]+=1.0;
      }
    }
    else {
      if (myrank!=0) {
    	exchange_data[0]=-1;
    	exchange_data[1]=-1;
    	exchange_data[2]=-1;
    	exchange_data[3]=-1;
      }
    
      if (myrank<your_rank)
    	printf("%d-th Rej %d-%d : /_\\ = %8.4lf \n EAAmi=%8.4lf ECGmi=%8.4lf EZmi=%8.4lf EAAni=%8.4lf ECGni=%8.4lf EZni=%8.4lf\n EAAmj=%8.4lf ECGmj=%8.4lf EZmj=%8.4lf EAAnj=%8.4lf ECGnj=%8.4lf EZnj=%8.4lf\n KZAAm=%8.4lf KZAAn=%8.4lf KZCGm=%8.4lf KZCGn=%8.4lf\n",
    	       i+1,myrank,your_rank,myrank,delta,
    	       EAAm_Xi,ECGm_Xi,EZm_Xi,EAAn_Xi,ECGn_Xi,EZn_Xi,
    	       EAAm_Xj,ECGm_Xj,EZm_Xj,EAAn_Xj,ECGn_Xj,EZn_Xj,
    	       ZData.KZAA,ZData.KZCG,KZAAn,KZCGn);
    }
    
    if (myrank!=0) MPI_Send(&exchange_data, 4, MPI_INT, 0, tag, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (myrank==0) {
      for (j=1;j<numRE;++j) {
    	MPI_Recv(&exchange_data_receiv, 4, MPI_INT, j, tag, MPI_COMM_WORLD, status);
    	exchange_data_sum[j*4+0]=exchange_data_receiv[0];
    	exchange_data_sum[j*4+1]=exchange_data_receiv[1];
    	exchange_data_sum[j*4+2]=exchange_data_receiv[2];
    	exchange_data_sum[j*4+3]=exchange_data_receiv[3];
    	m=exchange_data_receiv[2];
    	n=exchange_data_receiv[3];
    	//	printf("%d:210 j=%d m=%d n=%d\n",myrank,j,m,n);
    	if ( m != -1 && n != -1 ) acc_ratio[m][n]+=1.0;
      }
      for (j=1;j<numRE;++j) {
    	REMD_exchange_purmutation_funcs(exchange_data_sum[j*4+0],exchange_data_sum[j*4+1],
    					exchange_data_sum[j*4+2],exchange_data_sum[j*4+3]);
      }
    
      for (j=0;j<numRE;++j) exchanged_data[j]=index_parameters[j];
      for (j=0;j<numRE;++j) exchanged_data[j+numRE]=index_replicas[j];
    
      for (j=1;j<numRE;++j) {
    	MPI_Send(exchanged_data, numRE*2, MPI_INT, j, tag, MPI_COMM_WORLD);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (myrank!=0) {
      MPI_Recv(exchanged_data_receiv, numRE*2, MPI_INT, 0, tag, MPI_COMM_WORLD, status);
      for (j=0;j<numRE;++j) {
    	index_parameters[j]=exchanged_data_receiv[j];
    	index_replicas[j]=exchanged_data_receiv[numRE+j];
      }
    }
    
    /******************************************************************************************************/
    /* if (myrank==1) {											  */
    /*   printf("%d: 244 \n",myrank);									  */
    /*   for (j=0;j<numRE;++j) {									  */
    /* 	printf("i=%d-m=%d i=%d-m=%d ",j,REMD_purmutation_func(j),j,index_parameters[j]);		  */
    /* 	printf("m=%d-i=%d m=%d-i=%d\n",j,REMD_purmutation_inverse_func(j),j,index_replicas[j]);		  */
    /*   }												  */
    /*   printf("%d: 249 \n",myrank);									  */
    /* }												  */
    /******************************************************************************************************/
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  for (i=0;i<numRE;++i) for (j=0;j<numRE;++j) acc_ratio[i][j]=acc_ratio[i][j]/numEX*2.0;
  
  return acc_ratio;
}
