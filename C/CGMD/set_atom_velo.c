#include <stdio.h>
#include <math.h>

#include "ABA.h" // #include "ABA_multimer.h"
#include "gener.h"
#include "physics.h"
#include "MD.h"

// 各原子の速度の計算を行う関数
void set_atom_velo(void)
{
	int i;
	int alpha;

	for(i = 0; i < prot.num_atom; ++i)
	{
		for(alpha = 0; alpha < 3; ++alpha)
		{
		  //			prot.velo[i][alpha]/*A/s*/ 
		  //	= (  prot.coord[i][alpha]/*A*/
		  //	   - prot.old_coord[i][alpha]/*A*/)/deltat/**10e-10*/*1.0e-10*1.0e12;
		  //	prot.acc[i][alpha]/*A/s*/
		  //	=  (  prot.velo[i][alpha]/*A*/
		  //	    - prot.old_velo[i][alpha]/*A*/)/deltat/**10e-10*/*1.0e-10*1.0e12;
		}
	}
}
