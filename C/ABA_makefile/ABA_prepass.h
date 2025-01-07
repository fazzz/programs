#ifndef INCLUDE_ABA_P
#define INCLUDE_ABA_P

void ABAp_prepass(CLT *clt,double *qvel,int numclt,int numatom,double *crd);

void ABAp_prepass_TermOn(CLT *clt,double *qvel,int numclt,int numatom,double *crd,double *vel_Term);

void ABAp_Cofr(double Cofrc[6],double Spvel[6],double IM[6][6],double qCOM[3],double Sum_Mass);
void ABA_Coacc(double Coacc[6],double Spvel[6],double Spvel_P[6],double qvel,
	       int origin_atom_a,int origin_atom_a_Parent,int nNumAtomOfminousone,
	       double *crd,double trans_A_to_CN_P[3][3],double TMat[6][6]);
void ABAp_Spvelo(double Spvel[6],double qvel,double TMat[6][6],double SpvelPart[6]);

#endif

