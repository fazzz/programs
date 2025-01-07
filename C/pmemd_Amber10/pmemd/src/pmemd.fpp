#include "copyright.i"

!*******************************************************************************
!
! Program:      PMEMD 10.0, based on sander version 10
!
! Description:  The Molecular Dynamics/NMR Refinement/Modeling Module of the
!               AMBER Package, high performance version.
! Author:       PMEMD 10.0 has been developed entirely by Robert E. Duke of the
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

  implicit none

! Local variables:

  double precision      :: max_erfc_relerr
  integer               :: new_limit    ! new stack limit
  integer               :: i
  integer               :: num_ints = 0
  integer               :: num_reals = 0

  call second(run_start_cputime)
  call wall(run_start_walltime)

! Set up parallel execution environment:

#ifdef MPI
  call mstartup(mytaskid, numtasks)

  master = mytaskid .eq. 0 ! Make task 0 the master:
  if (numtasks .lt. 2 .and. master) then
    write(mdout, *) &
      'MPI version of PMEMD must be used with 2 or more processors!'
    call mexit(6, 1)
  end if
#else
  master = .true.  ! In the single-threaded version, the 1 process is master
#endif

! Reset the stack limits if you can:

  call unlimit_stack(new_limit)

  if (master .and. new_limit .gt. 0) then
    write(mdout, '(a,a,i10,a)') warn_hdr, &
      'Stack usage limited by a hard resource limit of ', new_limit, ' bytes!'
    write(mdout, '(a,a)') extra_line_hdr, &
      'If segment violations occur, get your sysadmin to increase the limit.'
  end if

#ifdef TIME_TEST
  call init_test_timers   ! For mpi performance monitoring
  call enable_test_timers
#endif

! Uniprocessor and mpi master initialization.  The master does all i/o.

  if (master) call master_setup(num_ints, num_reals)

#ifdef MPI
  ! Do generic broadcasts of data used in pretty much all circumstances...
  call bcast_mdin_ctrl_dat
  call bcast_amber_prmtop_dat
  call bcast_inpcrd_dat(natom)
  call bcast_constraints_dat(natom, ibelly, ntr)
  call bcast_extra_pnts_nb14_dat
  call parallel_dat_setup(natom, num_ints, num_reals)
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
      write(mdout, '(a,i8,/)') '| Nonbonded Pairs Initial Allocation:', &
                               ipairs_maxsize
    end if
  end if

#ifdef MPI
  if (master) then
    write(mdout, '(a,i4,a,/)') '| Running AMBER/MPI version on ', numtasks, &
                               ' nodes'
    write(mdout, '(a)') ' '
  end if
#endif

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
! In this implementation of pmemd, when running molecular dynamics,
! parallelization occurs at the level of subroutine runmd(), and when running
! minimizations, parallelization occurs at the level of subroutine force().

  if (nmropt .ne. 0) call bcast_nmr_dat

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
      call mexit(0, 0)
    endif
  end if
#else
  call second(run_setup_end_cputime)
  call wall(run_setup_end_walltime)
#endif /* end MPI */

! Now do the dynamics or minimization:

  if (imin .eq. 0) then         ! Do molecular dynamics:

#ifdef MPI
    call runmd(natom, atm_crd, atm_mass, atm_frc, atm_vel, atm_last_vel, &
               gbl_my_atm_lst)
#else
    call runmd(natom, atm_crd, atm_mass, atm_frc, atm_vel, atm_last_vel)
#endif /* MPI */
  else                          ! Do minimization:
    call runmin_master(natom, atm_crd, atm_frc, atm_igraph)
    ! Write restart file. atm_vel is not actually used.
    if (master) call write_restart(16, natom, atm_crd, atm_vel, 0.d0)
  endif

#ifdef MPI
! If doing minimization, set and broadcast notdone to inform other nodes that
! we are finished calling force()

  if (imin .ne. 0) then
    notdone = 0
    call mpi_bcast(notdone, 1, mpi_integer, 0, mpi_comm_world, err_code_mpi)
  endif

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

#ifdef TIME_TEST
  call print_test_timers        ! Debugging output of performance timings
#endif
  call mexit(6, 0)

end program pmemd
