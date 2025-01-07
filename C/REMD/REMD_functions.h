
#ifndef INCLUDE_REMD_fs
#define INCLUDE_REMD_fs

int *index_replicas;
int *index_parameters;

int REMD_purmutation_func(int i);
int REMD_purmutation_inverse_func(int m);
int REMD_exchange_purmutation_funcs(int i,int j,int m, int n);
int REMD_ini_purmutation_funcs(int num);

#endif
