#define _GNU_SOURCE  

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <netcdf.h>
#include <getopt.h>

#include "TOPO.h"
#include "LA.h"
#include "mymath.h"

#include "PTL.h"
#include "EF.h"

#define NVT 1
#define NVE 0

#define ON 1
#define OFF 0

void usage(char *progname);

int main(int argc, char *argv[]) {
  int i,j,k,l;
  int ii,jj,kk,ll;

  double fi[3],fj[3],fk[3],fl[3];
  double m[3],n[3],m_n[3],n_n[3],lm,ln;
  double vij[3],vkj[3],vkl[3];
  double lkj;
  double vijvkj,vklvkj;

  double dvdpsi;
  
  double atom[4][3];
  double dihed;
  double p_d_t=0.0;

  double pi;

  pi=acos(-1.0);

  atom[0][0]=16.01212;    atom[0][1]=1.78230;   atom[0][2]=-8.98534;
  atom[1][0]=16.46528;    atom[1][1]=3.11455;   atom[1][2]=-8.94461;
  atom[2][0]=17.22421;    atom[2][1]=3.53744;   atom[2][2]=-7.70562;
  atom[3][0]=17.49137;    atom[3][1]=2.72338;   atom[3][2]=-6.64599;
  
  for (j=0;j<3;++j) {
    vij[j] = atom[1][j]-atom[0][j];
    vkj[j] = atom[1][j]-atom[2][j];
    vkl[j] = atom[3][j]-atom[2][j];
  }
  lkj=sqrt(inprod(vkj,vkj,3));
      
  outprod(vij,vkj,m);
  outprod(vkj,vkl,n);
  lm=sqrt(inprod(m,m,3));
  ln=sqrt(inprod(n,n,3));
  for (j=0;j<3;++j) {
    m_n[j]=m[j]/lm;
    n_n[j]=n[j]/ln;
  }

  dihed=acos(inprod(m_n,n_n,3));
  if (inprod(vij,n,3)>0) dihed=-dihed;
  if (dihed<0.0) dihed=2.0*pi+dihed;
  
  vijvkj=inprod(vij,vkj,3);
  vklvkj=inprod(vkl,vkj,3);
  
  //  dvdpsi=Kd1*sin(dihed-dih_equ[i])+3.0*Kd2*sin(3.0*(dihed-dih_equ[i]));
  
  //  p_d[i] = Kd1*(1.0-cos(dihed-dih_equ[i]))+Kd2*(1.0-cos(3.0*(dihed-dih_equ[i])));
  //  p_d_t += p_d[i];
  
  /*************************************************************************/
  /* for (j=0;j<3;++j) {						   */
  /*   fi[j] = -dvdpsi*lkj*m[j]/(lm*lm);				   */
  /*   fl[j] =  dvdpsi*lkj*n[j]/(ln*ln);				   */
  /*   fj[j] = (-fi[j]+(vijvkj/(lkj*lkj))*fi[j]-(vklvkj/(lkj*lkj))*fl[j]); */
  /*   fk[j] = (-fl[j]-(vijvkj/(lkj*lkj))*fi[j]+(vklvkj/(lkj*lkj))*fl[j]); */
  /*   									   */
  /*   f_d[ii][j] += fi[j]*4.184070*100.0;				   */
  /*   f_d[jj][j] += fj[j]*4.184070*100.0;				   */
  /*   f_d[kk][j] += fk[j]*4.184070*100.0;				   */
  /*   f_d[ll][j] += fl[j]*4.184070*100.0;				   */
  /* }									   */
  /*************************************************************************/

  return 0.0;
}
