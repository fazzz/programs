/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /home/case/cvsroot/amber11/AmberTools/src/netcdf/src/ncgen/init.c,v 9.1 2007/11/15 21:44:48 jmongan Exp $
 *********************************************************************/

#include <stdio.h>
#include <netcdf.h>
#include "generic.h"
#include "ncgen.h"
#include "genlib.h"

extern int netcdf_flag;
extern int c_flag;
extern int fortran_flag;

struct dims *dims;		/* table of netcdf dimensions */

int ncid;			/* handle for netCDF */
int ndims;			/* number of dimensions declared for netcdf */
int nvars;			/* number of variables declared for netcdf */
int natts;			/* number of attributes */
int nvdims;			/* number of dimensions for variables */
int dimnum;			/* dimension number index for variables */
int varnum;			/* variable number index for attributes */
int valnum;			/* value number index for attributes */
int rec_dim;			/* number of the unlimited dimension, if any */
size_t var_len;			/* variable length (product of dimensions) */
size_t rec_len;			/* number of elements for a record of data */
size_t var_size;		/* size of each element of variable */

struct vars *vars;		/* a malloc'ed list */

struct atts *atts;		/* table of variable and global attributes */

void
init_netcdf(void) {			/* initialize global counts, flags */
    
    clearout();			/* reset symbol table to empty */
    ndims = 0;
    nvars = 0;
    rec_dim = -1;		/* means no unlimited dimension (yet) */
}
