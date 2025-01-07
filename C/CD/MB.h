
#define phi 1
#define psi 2
#define omega 3
#define kai1 4

int countdihed(int flagKOP);

double *CD(double *crd,int flagKOP);

double pick_dihed(  double atom_i[3],double atom_j[3],double  atom_k[3],double atom_l[3], int flag, double old_value);

double pick_bond_leng(  double atom_i[3],double atom_j[3]);

int check_phi_psi_omega_kai1(int atom1,int atom2, int atom3, int atom4 );

int pick_atom_dihed_pairs(int *atom_dihed_pair,int flag);

void pick_dihed_all(double *coord,double *dihed,int numdihed,int *atom_dihed_pair);

double pick_angle(double atom_i[3],double atom_j[3],double atom_k[3],int flag,double old_value);

double wraped_angle(double angle_rad, double pi);
