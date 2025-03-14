/*********************************************************************
 * Copyright 1993-2006, UCAR/Unidata
 * See COPYRIGHT file for copying and redistribution conditions.
 *
 * Test driver for netCDF implementation.  This program performs tests
 * against the netCDF specification for all user-level functions in an
 * implementation of the netCDF library.  Must be invoked from a
 * directory in which the invoker has write permission.
 *
 * $Id: driver.c,v 9.1 2007/11/15 21:44:49 jmongan Exp $
 *********************************************************************/

#include <config.h>
#include <netcdf.h>
#include <stdio.h>
#include "tests.h"

/* Test everything for classic and 64-bit offsetfiles. If netcdf-4 is
 * included, that means another whole round of testing. */
#ifdef USE_NETCDF4
#define NUM_FORMATS (3)
#else
#define NUM_FORMATS (2)
#endif

int
main()
{
    extern int ncopts;		/* netCDF error options */
    static char *testfiles[] = {"nonesuch", "nctest_classic.nc", 
				"nctest_64bit_offset.nc", "nctest_netcdf4.nc"};
    char *testfile;
    int i, nerrs = 0;

    ncopts &= ~NC_FATAL;	/* make errors nonfatal */
    ncopts &= ~NC_VERBOSE;	/* turn off error messages */

    fprintf(stderr, "Testing V2 API with %d different netCDF formats.\n", 
	   NUM_FORMATS);

    for (i = 1; i <= NUM_FORMATS; i++)
    {
       switch (i) 
       {
	  case NC_FORMAT_CLASSIC:
	     nc_set_default_format(NC_FORMAT_CLASSIC, NULL);
	     fprintf(stderr, "\n\nSwitching to netCDF classic format.\n");
	     break;
	  case NC_FORMAT_64BIT:
	     nc_set_default_format(NC_FORMAT_64BIT, NULL);
	     fprintf(stderr, "\n\nSwitching to 64-bit offset format.\n");
	     break;
#ifdef USE_NETCDF4 
	  case NC_FORMAT_NETCDF4: /* actually it's _CLASSIC. */
	     nc_set_default_format(NC_FORMAT_NETCDF4_CLASSIC, NULL);
	     fprintf(stderr, "\n\nSwitching to netCDF-4 format (with NC_STRICT_NC3).\n");
	     break;
#endif
	  default:
	     fprintf(stderr, "Unexpected format!\n");
	     return 2;
       }
       testfile = testfiles[i];

       /* Run all the tests for this format. */
       nerrs += test_nccreate(testfile);
       nerrs += test_ncopen(testfile);
       nerrs += test_ncredef(testfile);
       nerrs += test_ncendef(testfile);
       nerrs += test_ncclose(testfile);
       nerrs += test_ncinquire(testfile);
       nerrs += test_ncsync(testfile);
       nerrs += test_ncabort(testfile);
       nerrs += test_ncdimdef(testfile);
       nerrs += test_ncdimid(testfile);
       nerrs += test_ncdiminq(testfile);
       nerrs += test_ncdimrename(testfile);
       nerrs += test_ncvardef(testfile);
       nerrs += test_ncvarid(testfile);
       nerrs += test_ncvarinq(testfile);
       nerrs += test_ncvarput1(testfile);
       nerrs += test_ncvarget1(testfile);
       nerrs += test_ncvarput(testfile);
       nerrs += test_ncvarget(testfile);
       nerrs += test_ncvarputg(testfile);
       nerrs += test_ncvargetg(testfile);
       nerrs += test_ncrecinq(testfile);
       nerrs += test_ncrecput(testfile);
       nerrs += test_ncrecget(testfile);
       nerrs += test_ncvarrename(testfile);
       nerrs += test_ncattput(testfile);
       nerrs += test_ncattinq(testfile);
       nerrs += test_ncattget(testfile);
       nerrs += test_ncattcopy(testfile, "test2.nc");
       /* These tests don't work with HDF5 1.8.0 alpha releases (yet)
	* but are expected to one day. The --enable-hdf5-alpha-release
	* option to configure allows the user to turn off these tests
	* for netCDF-4 files. */
#ifdef HDF5_ALPHA_RELEASE
       if (i != NC_FORMAT_NETCDF4)
#endif	  
       {
	  nerrs += test_ncattname(testfile);
	  nerrs += test_ncattrename(testfile);
	  nerrs += test_ncattdel(testfile);
       }
       nerrs += test_nctypelen();
    }
    
    fprintf(stderr, "\nTotal number of failures: %d\n", nerrs);

    if (nerrs)
    {
       fprintf(stderr, "nctest FAILURE!!!\n");
       return 2;
    }
    else
       fprintf(stderr, "nctest SUCCESS!!!\n");
    return 0;
}
