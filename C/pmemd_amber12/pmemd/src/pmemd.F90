#include "copyright.i"

!*******************************************************************************
!
! Program:      PMEMD 11.0, originally based on sander version 6
!
! Description:  The Molecular Dynamics/NMR Refinement/Modeling Module of the
!               AMBER Package, high performance version.
! Author:       PMEMD 11.0 has been developed by Robert E. Duke of the
!               University of North Carolina-Chapel Hill Chemistry Department
!               and NIEHS.  The development was originally based on Amber 6.0
!               sander source code.  The product is intended to be a fast
!               implementation of the Particle Mesh Ewald method for molecular
!               dynamics and minimizations.  The code now also supports
!               generalized Born and ALPB molecular dynamics and minimizations.
!               This work was done with the support of Prof. Lee Pedersen and
!               his group at the University of North Carolina-Chapel Hill,
!               Tom Darden at NIEHS, and the Amber core development group.
!               Funding support was provided by NIH grant HL-06350 (PPG),
!               NSF grant 2001-0759-02 (ITR/AP), and intramural NIH funds.
!
!               Extensions to support the CHARMM force field were made by
!               Mike Crowley, Mark Williamson, and Ross Walker.
!
!               NVIDIA GPU acceleration support was written by Scott Le Grand,
!               Duncan Poole and Ross Walker.
!
!               Replica exchange support, igb=8, and GB/SA written by 
!               Jason Swails, Adrian Roitberg, and Ross Walker, 
!               Funded under NSF SSE Grant to RCW & AER
!*******************************************************************************

program pmemd

  use gb_alltasks_setup_mod
  use pme_alltasks_setup_mod
  use constraints_mod
  use extra_pnts_nb14_mod
#ifdef DIRFRC_EFS
  use ene_frc_splines_mod
#endif /* DIRFRC_EFS */
  use gb_ene_mod
  use gbl_constants_mod
  use gbl_datatypes_mod
  use pme_setup_mod
  use img_mod
  use inpcrd_dat_mod
  use master_setup_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
#ifdef MPI
  use multipmemd_mod, only : setup_groups, free_comms
#endif
  use nb_pairlist_mod
  use nmr_calls_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use random_mod
  use runfiles_mod
  use runmd_mod
  use runmin_mod
  use shake_mod
  use timers_mod
  use bintraj_mod
  use pbc_mod
  use nbips_mod
#ifdef MPI
  use pme_force_mod,  only : alloc_force_mem, dealloc_force_mem
  use remd_mod,       only : remd_method, bcast_remd_method, remd_setup, &
                             slave_remd_setup, remd_cleanup
  use remd_exchg_mod, only : setup_remd_randgen
#endif
  use amd_mod
  use file_io_mod,       only : amopen

  implicit none

! Local variables:

  double precision      :: max_erfc_relerr
  integer               :: new_stack_limit    ! new stack limit
  integer               :: i
  integer               :: num_ints = 0
  integer               :: num_reals = 0
  logical               :: terminal_flag = .false.

  call second(run_start_cputime)
  call wall(run_start_walltime)

#ifdef MPI
! Establish pmemd communicators and process ranks
  call setup_groups

  master = mytaskid .eq. 0 ! Make task 0 the master:

! CUDA version with GB can run on 1 GPU / task.
! All other combinations require .gt. 1
#ifndef CUDA
  if (numtasks .lt. 2 .and. master) then
    write(mdout, *) &
      'MPI version of PMEMD must be used with 2 or more processors!'
    call mexit(6, 1)
  end if
#endif

#else
  master = .true.  ! In the single-threaded version, the 1 process is master
#endif

! Reset the stack limits if you can:

  call unlimit_stack(new_stack_limit)

#ifdef TIME_TEST
  call init_test_timers   ! For mpi performance monitoring
  call enable_test_timers
#endif

! Create gpu context
#ifdef CUDA
#ifdef MPI
  call gpu_startup(mytaskid, numtasks, pmemd_comm_number)
#else
  call gpu_startup()
#endif
#endif

  if (master) then
    call master_setup(num_ints, num_reals, new_stack_limit, terminal_flag)
#ifdef CUDA
  else
    !RCW: Call to set device is mostly redundant now due to
    !     removal of -gpu command line argument. But leave for now.
    call gpu_set_device(-1)
    call gpu_init()
#ifdef MPI
    call gpu_send_slave_device_info()
#endif    
#endif
  end if

#ifdef MPI
  call bcast_logical(terminal_flag, pmemd_comm)
  if (terminal_flag) call mexit(6,0)
  ! Do generic broadcasts of data used in pretty much all circumstances...
  call bcast_remd_method
  call bcast_mdin_ctrl_dat
  call bcast_amber_prmtop_dat
  call bcast_inpcrd_dat(natom)
  call bcast_constraints_dat(natom, ibelly, ntr)
  call bcast_extra_pnts_nb14_dat
  call parallel_dat_setup(natom, num_ints, num_reals)
  call alloc_force_mem(natom,num_reals,ips)
#ifdef CUDA
  if (ntb .ne. 0) then
    call bcast_pbc
    call bcast_mdin_ewald_dat
    if (.not. master) then
      call gpu_init_pbc(pbc_box(1), pbc_box(2), pbc_box(3), pbc_alpha, &
                        pbc_beta, pbc_gamma, uc_volume, uc_sphere, &
                        vdw_cutoff + skinnb, pbc_box, reclng, cut_factor, &
                        ucell, recip)
    end if
  end if
  if (nmropt .ne. 0) call bcast_nmr_dat  
#endif  /* CUDA */
  ! Set up REMD if we're actually performing REMD
  if (remd_method .ne. 0) then
    if (master) then
      call remd_setup(numexchg)
    else
      call slave_remd_setup
    end if
    call setup_remd_randgen
  end if
#endif /* MPI */

  ! Set up AMD 
  if (iamd .gt. 0) then
    call amd_setup(ntwx)
  endif

! If a terminal flag was put on the command-line (--help, --version, etc.), then
! just exit here
if (terminal_flag) call mexit(6,0)

! Initialize system
#ifdef CUDA
  call gpu_setup_system(natom, tol, ntf, ntb, ips, ntp, ntt, vrand)
  call gpu_upload_crd(atm_crd)
  call gpu_upload_charges(atm_qterm)
  call gpu_upload_masses(atm_mass)
  call gpu_upload_frc(atm_frc)
  call gpu_upload_vel(atm_vel)
  call gpu_upload_last_vel(atm_last_vel)
  call gpu_init_extra_pnts_nb14(gbl_frame_cnt, ep_frames, ep_lcl_crd)
  call gpu_constraints_setup(natc, atm_jrc, atm_weight, atm_xc)
#endif

! The following call does more uniprocessor setup and mpi master/slave setup.

  if (using_pme_potential) then
    call pme_alltasks_setup(num_ints, num_reals)
  else if (using_gb_potential) then
    call gb_alltasks_setup(num_ints, num_reals)
  end if



! Call shake_setup, which will tag those bonds which are part of 3-point
! water molecules and also set up data structures for non-fastwater shake.
! Constraints will be done for waters using a fast analytic routine -- dap.
! Currently, ALL processes call shake_setup, as it is now very fast, and this
! is probably cheaper than broadcasting the results (also, there is now atom
! selection in shake setup under mpi).

  call shake_setup(num_ints, num_reals)

  if (using_pme_potential) then
    call final_pme_setup(num_ints, num_reals)
  else if (using_gb_potential) then
    call final_gb_setup(natom, num_ints, num_reals)
  end if

! Deallocate data that is not needed after setup:

  if (using_pme_potential) then
    num_ints = num_ints - size(atm_numex) - size(gbl_natex)
    deallocate(atm_numex, gbl_natex)
  end if

  if (imin .ne. 0) then
    num_reals = num_reals - size(atm_vel)
    deallocate(atm_vel)
  end if
  
! Should be all done with setup. Report dynamic memory allocation totals:

  if (master) then
    write(mdout, '(a)')       '| Dynamic Memory, Types Used:'
    write(mdout, '(a,i14)')   '| Reals      ', num_reals
    write(mdout, '(a,i14,/)') '| Integers   ', num_ints
    if (using_pme_potential) then
      write(mdout, '(a,i12,/)') '| Nonbonded Pairs Initial Allocation:', &
                               ipairs_maxsize
    end if
  end if

#ifdef CUDA
#ifdef MPI
    if (master) then
#endif
      call gpu_write_memory_info(mdout)
#ifdef MPI
    end if
#endif    
#endif

#ifdef MPI
  if (master) then
    write(mdout, '(a,i4,a,/)') '| Running AMBER/MPI version on ', numtasks, &
                               ' nodes'
    write(mdout, '(a)') ' '
  end if
#endif

! Prepare for Isotropic periodic sum of nonbonded interaction
   if (ips .gt. 0)then
     call ipssys(atm_crd)
   endif

  if (master) write(mdout,'(80(''-'')/''   4.  RESULTS'',/80(''-'')/)')

  if (using_pme_potential) then
    call chk_switch(max_erfc_relerr)
#ifdef DIRFRC_EFS
    call chk_ef_spline_tables(max_erfc_relerr)
#endif /* DIRFRC_EFS */
  end if

! We reset the random number generator here to be consistent with sander 10.
! I think this is pretty much unnecessary in pmemd, but we will only get
! sander-consistent results if we do it.

  call amrset(ig + 1)


#ifdef MPI
! If we're running replica exchange, we need to reset our random number
! generator to match the random streams that sander produces so we can
! reproduce those results here

  if (remd_method .ne. 0) call amrset(ig + 17 * (repnum-1))

! In this implementation of pmemd, when running molecular dynamics,
! parallelization occurs at the level of subroutine runmd(), and when running
! minimizations, parallelization occurs at the level of subroutine force().

#ifndef CUDA
  if (nmropt .ne. 0) call bcast_nmr_dat
#endif

  call second(run_setup_end_cputime)
  call wall(run_setup_end_walltime)

  ! Parallelization of minimization, nonmaster nodes:

  if (imin .ne. 0) then
    if (.not. master) then         ! All nodes do only force calc
      call runmin_slave(natom, atm_crd, atm_frc, gbl_img_atm_map)
      call second(run_end_cputime)
      call profile_cpu(imin, igb, loadbal_verbose)

#ifdef TIME_TEST
      call print_test_timers    ! Debugging output of performance timings
#endif
      call free_comms
      call mexit(0, 0)
    endif
  end if
#else
  call second(run_setup_end_cputime)
  call wall(run_setup_end_walltime)
#endif /* end MPI */

  ! Initialize the printing of ongoing time and performance summaries. We call this
  ! here after all the setup is done so we don't end up including all the startup
  ! time etc.
  if (master) call print_ongoing_time_summary(0,0,0.0d0,0)

! Now do the dynamics or minimization:

  if (imin .eq. 0) then         ! Do molecular dynamics:

#ifdef MPI
    call runmd(natom, atm_crd, atm_mass, atm_frc, atm_vel, atm_last_vel, &
               gbl_my_atm_lst, remd_method, numexchg)
#else
    call runmd(natom, atm_crd, atm_mass, atm_frc, atm_vel, atm_last_vel)
#endif /* MPI */
  else                          ! Do minimization:
    call runmin_master(natom, atm_crd, atm_frc, atm_igraph)
    ! Write restart file. atm_vel is not actually used.
    call write_restart(restrt, natom, atm_crd, atm_vel, &
                       0.d0, .true., restrt_name)
    close (restrt)
  endif

#ifdef MPI
! If doing minimization, set and broadcast notdone to inform other nodes that
! we are finished calling force()

  if (imin .ne. 0) then
    notdone = 0
    call mpi_bcast(notdone, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)
  endif

  call dealloc_force_mem(ips)

  if (imin .eq. 0) then
    call remd_cleanup
  end if

  if (master) then 
    if (ioutfm .eq. 1) then
      call close_binary_files
    else
      if (ntwx .gt. 0) close(mdcrd)
      if (ntwv .gt. 0) close(mdvel)
    end if
    if (ntwe .gt. 0) close(mden)
  else

    call second(run_end_cputime)
    call profile_cpu(imin, igb, loadbal_verbose)

#ifdef TIME_TEST
    call print_test_timers      ! Debugging output of performance timings
#endif
#ifdef CUDA
! Shut down GPU
  call gpu_shutdown()
#endif
    call free_comms
    call mexit(0, 0)
  endif
#endif /* end MPI */

  ! Master prints timings (averages for mpi) in mdout:

  call second(run_end_cputime)
  call wall(run_end_walltime)

  write(mdout,'(80(1H-)/,''   5.  TIMINGS'',/80(1H-)/)')
#ifdef MPI
  call profile_cpu(imin, igb, loadbal_verbose)
#else
  call profile_cpu(imin, igb, 0)
#endif

  !Write final performance numbers to mdout
  call print_ongoing_time_summary(nstlim,nstlim,dt,mdout)

#ifdef MPI
  write(mdout, '(/, a, f11.2, a)') &
        '|  Master Setup CPU time:     ', &
        run_setup_end_cputime - run_start_cputime, ' seconds'
  write(mdout, '(   a, f11.2, a)') &
        '|  Master NonSetup CPU time:  ', &
        run_end_cputime - run_setup_end_cputime, ' seconds'
  write(mdout, '(   a, f11.2, a, f9.2, a)') &
        '|  Master Total CPU time:     ', &
        run_end_cputime - run_start_cputime, ' seconds', &
        (run_end_cputime - run_start_cputime) / 3600.d0, ' hours'
  write(mdout, '(/, a, i9, a)') &
        '|  Master Setup wall time:   ', &
        run_setup_end_walltime - run_start_walltime, '    seconds'
  write(mdout, '(   a, i9, a)') &
        '|  Master NonSetup wall time:', &
        run_end_walltime - run_setup_end_walltime, '    seconds'
  write(mdout, '(   a, i9, a, f9.2, a)') &
        '|  Master Total wall time:   ', &
        run_end_walltime - run_start_walltime, '    seconds', &
        dble(run_end_walltime - run_start_walltime) / 3600.d0, ' hours'
#else
  write(mdout, '(/, a, f11.2, a)') &
        '|  Setup CPU time:     ', &
        run_setup_end_cputime - run_start_cputime, ' seconds'
  write(mdout, '(   a, f11.2, a)') &
        '|  NonSetup CPU time:  ', &
        run_end_cputime - run_setup_end_cputime, ' seconds'
  write(mdout, '(   a, f11.2, a, f9.2, a)') &
        '|  Total CPU time:     ', &
        run_end_cputime - run_start_cputime, ' seconds', &
        (run_end_cputime - run_start_cputime) / 3600.d0, ' hours'
  write(mdout, '(/, a, i9, a)') &
        '|  Setup wall time:   ', &
        run_setup_end_walltime - run_start_walltime, '    seconds'
  write(mdout, '(   a, i9, a)') &
        '|  NonSetup wall time:', &
        run_end_walltime - run_setup_end_walltime, '    seconds'
  write(mdout, '(   a, i9, a, f9.2, a)') &
        '|  Total wall time:   ', &
        run_end_walltime - run_start_walltime, '    seconds', &
        dble(run_end_walltime - run_start_walltime) / 3600.d0, ' hours'
#endif

  close(mdout)

#ifdef CUDA
! Shut down GPU
  call gpu_shutdown()
#endif



#ifdef TIME_TEST
  call print_test_timers        ! Debugging output of performance timings
#endif
#ifdef MPI
  call free_comms               ! free MPI communicators
#endif
  call mexit(6, 0)

end program pmemd
