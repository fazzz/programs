#ifndef INCLUDE_ABAb_P
#define INCLUDE_ABAb_P

void ABAbp_prepass(CLTb *clt,double *qvel,int numclt,int numatom,double *crd);

void ABAbp_prepass_TermOn(CLTb *clt,double *qvel,int numclt,int numatom,double *crd,double *vel_Term);

void ABAbp_Cofr(double Cofrc[6],double Spvel[6],double IM[6][6],double qCOM[3],double Sum_Mass);
void ABAb_Coacc(double Coacc[6],double Spvel[6],double Spvel_P[6],double qvel,
	       int origin_atom_a,int origin_atom_a_Parent,int nNumAtomOfminousone,
	       double *crd,double trans_A_to_CN_P[3][3],double TMat[6][6]);
void ABAbp_Spvelo(double Spvel[6],double qvel,double TMat[6][6],double SpvelPart[6]);

#endif

