#ifndef INCLUDE_MB
#define INCLUDE_MB

#define phi 1
#define psi 2
#define omega 3
#define kai1 4

double pick_dihed(  double atom_i[3],double atom_j[3],double  atom_k[3],double atom_l[3], 
		    int flag, double old_value);

double pick_angle(double atom_i[3],double atom_j[3],double atom_k[3],int flag,double old_value);

double pick_bond_leng(  double atom_i[3],double atom_j[3]);

#endif
