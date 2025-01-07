#
#************************************************************************
#                              AMBER                                   **
#                                                                      **
#               Copyright (c) 1986, 1991, 1995, 1997, 1999, 2004, 2008 **
#                Regents of the University of California               **
#                       All Rights Reserved.                           ** 
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************

# Makefile for Amber 12  PMEMD 
SHELL=/bin/sh
AMBER=Amber12
CONFIG_FILE=config.h
CONFIG_COMMAND=./configure
PARALLEL_TOKEN=DMPI
CUDA_TOKEN=DCUDA
CUDA_DPDP_TOKEN=Duse_DPDP
CUDA_SPSP_TOKEN=Duse_SPSP

# Platform-specific info should be found in ../../config.h
include $(CONFIG_FILE)

OBJS=   gbl_constants$(OSFX) gbl_datatypes$(OSFX) state_info$(OSFX) file_io_dat$(OSFX) \
        mdin_ctrl_dat$(OSFX) mdin_ewald_dat$(OSFX) mdin_debugf_dat$(OSFX) prmtop_dat$(OSFX) \
        inpcrd_dat$(OSFX) dynamics_dat$(OSFX) img$(OSFX) nbips$(OSFX) \
        parallel_dat$(OSFX) parallel$(OSFX) gb_parallel$(OSFX) \
        pme_direct$(OSFX) pme_recip_dat$(OSFX) pme_slab_recip$(OSFX) pme_blk_recip$(OSFX) \
        pme_slab_fft$(OSFX) pme_blk_fft$(OSFX) pme_fft_dat$(OSFX) fft1d$(OSFX) \
        bspline$(OSFX) pme_force$(OSFX) pbc$(OSFX) nb_pairlist$(OSFX) \
        nb_exclusions$(OSFX) cit$(OSFX) dynamics$(OSFX) bonds$(OSFX) angles$(OSFX) dihedrals$(OSFX) \
        extra_pnts_nb14$(OSFX) runmd$(OSFX) loadbal$(OSFX) shake$(OSFX) prfs$(OSFX) mol_list$(OSFX) \
        runmin$(OSFX) constraints$(OSFX) axis_optimize$(OSFX) gb_ene$(OSFX) veclib$(OSFX) gb_force$(OSFX) \
        timers$(OSFX) pmemd_lib$(OSFX) runfiles$(OSFX) file_io$(OSFX) bintraj$(OSFX)  pmemd_clib$(OSFX) \
        pmemd$(OSFX) random$(OSFX) degcnt$(OSFX) erfcfun$(OSFX) nmr_calls$(OSFX) nmr_lib$(OSFX) \
        get_cmdline$(OSFX) master_setup$(OSFX) pme_alltasks_setup$(OSFX) pme_setup$(OSFX) \
        ene_frc_splines$(OSFX) gb_alltasks_setup$(OSFX) nextprmtop_section$(OSFX) binrestart$(OSFX) \
        angles_ub$(OSFX) dihedrals_imp$(OSFX) cmap$(OSFX) charmm$(OSFX) charmm_gold$(OSFX) \
        findmask$(OSFX) remd$(OSFX) multipmemd$(OSFX) remd_exchg$(OSFX) amd$(OSFX)

install: configured_serial pmemd$(SFX)
	$(MV) pmemd$(SFX) .

parallel: configured_parallel pmemd.MPI$(SFX)
	$(MV) pmemd.MPI$(SFX) .

cuda: configured_cuda pmemd.cuda$(SFX)
	@( \
	  if grep $(CUDA_SPSP_TOKEN) $(CONFIG_FILE) 1> /dev/null 2>&1 ; then \
	    mv pmemd.cuda$(SFX) . ;\
	  elif grep $(CUDA_DPDP_TOKEN) $(CONFIG_FILE) 1> /dev/null 2>&1 ; then \
	    mv pmemd.cuda$(SFX) . ;\
	  else \
	    mv pmemd.cuda$(SFX) . ;\
	    ln -f -s pmemd.cuda_SPDP$(SFX) pmemd.cuda$(SFX);\
	  fi ;\
	)

cuda_parallel: configured_cuda configured_parallel pmemd.cuda.MPI$(SFX)
	@( \
	  if grep $(CUDA_SPSP_TOKEN) $(CONFIG_FILE) 1> /dev/null 2>&1 ; then \
	    mv pmemd.cuda.MPI$(SFX) pmemd.cuda_SPSP.MPI$(SFX) ;\
	  elif grep $(CUDA_DPDP_TOKEN) $(CONFIG_FILE) 1> /dev/null 2>&1 ; then \
	    mv pmemd.cuda.MPI$(SFX) pmemd.cuda_DPDP.MPI$(SFX) ;\
	  else \
	    mv pmemd.cuda.MPI$(SFX) pmemd.cuda_SPDP.MPI$(SFX) ;\
	    ln -f -s pmemd.cuda_SPDP.MPI$(SFX) pmemd.cuda.MPI$(SFX);\
	  fi ;\
	)

#Cannot use configured_win or configured_serial for Windows right now
#since this will not run in a dos shell.
#serial_win: configured_win configured_serial
#parallel_win: configured_win configured_parallel
install_win: pmemd$(SFX)
	$(MV) pmemd$(SFX) .

parallel_win: pmemd.MPI$(SFX)
	$(MV) pmemd.MPI$(SFX) .

pmemd$(SFX): $(NETCDF) $(OBJS)
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(LDOUT)$@ $(OBJS) -L$(LIBDIR) $(NETCDFLIB) $(LDFLAGS) $(PMEMD_FLIBSF)

pmemd.MPI$(SFX): $(NETCDF) $(OBJS)
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(LDOUT)$@ $(OBJS) -L$(LIBDIR) $(NETCDFLIB) $(LDFLAGS) $(PMEMD_FLIBSF)

pmemd.cuda$(SFX): $(NETCDF) $(OBJS) $(PMEMD_CU_LIBS)
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(PMEMD_CU_DEFINES) $(LDOUT)$@ $(OBJS) \
      $(PMEMD_CU_LIBS) -L$(LIBDIR) $(NETCDFLIB) $(LDFLAGS) $(PMEMD_FLIBSF)

pmemd.cuda.MPI$(SFX): $(NETCDF) $(OBJS) $(PMEMD_CU_LIBS)
	$(PMEMD_LD) $(PMEMD_FOPTFLAGS) $(PMEMD_CU_DEFINES) $(LDOUT)$@ $(OBJS) \
     $(PMEMD_CU_LIBS) -L$(LIBDIR) $(NETCDFLIB) $(LDFLAGS) $(PMEMD_FLIBSF)

$(INCDIR)/netcdf.mod:
	$(MAKE) -C ~/software/Amber12/amber12/src/netcdf/src install

depends:
	../makef90depends > f90depends

clean:
	$(RM) -f *.o *.obj *.mod pmemd.MPI$(SFX) pmemd$(SFX) pmemd.cuda$(SFX) pmemd.cuda.MPI$(SFX) *.d work.pc* 
	$(MAKE) -C ./cuda clean

# Control the suffixes available; this essentially eliminates built-in
# inference rules, done here to improve portability.
.SUFFIXES:
.SUFFIXES: .F90 .c $(OSFX)

.F90$(OSFX):
	$(PMEMD_F90) $(PMEMD_FOPTFLAGS) $(PMEMD_CU_DEFINES) $(NETCDFINC) -c $*.F90

.c$(OSFX):
	$(PMEMD_CC) $(PMEMD_COPTFLAGS) $(PMEMD_CU_DEFINES) $(NETCDFINC) -c $*.c

$(PMEMD_CU_LIBS):
	$(MAKE) -C ./cuda 

configured:
	@(if [ ! -f $(CONFIG_FILE) ] ; then \
	    echo "Error: $(CONFIG_COMMAND) must be executed before $(MAKE)!" ;\
	    exit 2 ;\ # $(CONFIG_COMMAND) ;\
	fi ;\
	)

configured_parallel: configured
	@(grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null || \
	    { echo "Error: $(CONFIG_FILE) is not of type parallel!" ;\
	      echo "  Rerun $(CONFIG_COMMAND) and specify an MPI implementation." ;\
	      exit 2 ;\
	    } ;\
	)

configured_serial: configured
	@(if grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
	    echo "Error: $(CONFIG_FILE) is of type parallel, not serial!" ;\
	    echo "  Rerun $(CONFIG_COMMAND) and specify -nopar." ;\
	    exit 2 ;\
	fi ;\
	)
	@(if grep $(CUDA_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
	    echo "Error: $(CONFIG_FILE) is of type cuda, not serial!" ;\
	    echo "  Rerun $(CONFIG_COMMAND) and specify -nopar." ;\
	    exit 2 ;\
	fi ;\
	)

configured_cuda: configured
	@(grep $(CUDA_TOKEN) $(CONFIG_FILE) > /dev/null || \
	    { echo "Error: $(CONFIG_FILE) is not of type cuda!" ;\
	      echo "  Rerun $(CONFIG_COMMAND) and specify a cuda implementation." ;\
	      exit 2 ;\
	    } ;\
	)

configured_win: configured
	@(grep $(WINDOWS_TOKEN) $(CONFIG_FILE) > /dev/null || \
	    { echo "Error: $(CONFIG_FILE) is not of type Windows" ;\
	      echo "  Rerun $(CONFIG_COMMAND) and specify -windows." ;\
	      exit 2 ;\
	    } ;\
	)

## Dependencies for f90 modules and include files generated by makef90depends:

include f90depends
