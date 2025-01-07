
#ifndef INCLUDE_PDBE
#define INCLUDE_PDB

#include "PDB.h"

int replace_ATOMNAME(PDBF PDB, char *ATOMNAME1, char *ATOMNAME2, int num);
int replaceall_ATOMNAME(PDBF PDB, char *ATOMNAME, int num);
int delete_ATOMNAME(PDBF PDB, char *ATOMNAME, int num);
int pickup_ATOMNAME(PDBF PDB, char *ATOMNAME,int num);
int put_occupancy(PDBF PDB,double occupancy,int num);
int show_line(PDBF PDB, int num);
int show_line_by_ATOMNAME(PDBF PDB,char *ATOMNAME, int num);
int show_line_by_RES(PDBF PDB,char *RESNAME, int num);


#endif
