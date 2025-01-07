#ifndef INCLUDE_TOPO
#define INCLUDE_TOPO

double len(double atom_i[3],double atom_j[3]);
double ang(double atom_i[3],double atom_j[3],double atom_k[3]);
double dih(double atom_i[3],double atom_j[3],double  atom_k[3],double atom_l[3]);
void csdih(double atom_i[3],double atom_j[3],double  atom_k[3],double atom_l[3], double *cs, double *sn);

#endif
