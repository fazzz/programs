#include "copyright.i"

!*******************************************************************************
!
! Module: runmd_mod
!
! Description: Sander 9-compatible MD.
!              
!*******************************************************************************

module runmd_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  runmd
!
! Description: <TBS>
!              
!*******************************************************************************

#ifdef MPI
subroutine runmd(atm_cnt, crd, mass, frc, vel, last_vel, my_atm_lst)
#else
subroutine runmd(atm_cnt, crd, mass, frc, vel, last_vel)
#endif

  use cit_mod
  use constraints_mod
  use degcnt_mod
  use dynamics_mod
  use dynamics_dat_mod
  use extra_pnts_nb14_mod
  use gb_force_mod
  use gb_ene_mod
  use pbc_mod
  use pme_force_mod
  use file_io_mod
  use file_io_dat_mod
  use img_mod
  use loadbal_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use nb_pairlist_mod
  use nmr_calls_mod
  use parallel_dat_mod
  use parallel_mod
  use gb_parallel_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use runfiles_mod
  use shake_mod
  use timers_mod
  use state_info_mod
  use bintraj_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: mass(atm_cnt)
  double precision      :: frc(3, atm_cnt)
  double precision      :: vel(3, atm_cnt)
  double precision      :: last_vel(3, atm_cnt)
#ifdef MPI
  integer               :: my_atm_lst(*)
#endif

! Local variables:

  double precision      :: aamass
  double precision      :: boltz2
  double precision      :: dtcp
  double precision      :: dtx, dtx_inv, half_dtx
  type(gb_pot_ene_rec)  :: gb_pot_ene
  type(pme_pot_ene_rec) :: pme_pot_ene
  double precision      :: pme_err_est
  double precision      :: eke
  double precision      :: ekmh
  double precision      :: ekpbs
  double precision      :: ekph
  double precision      :: ekin0
  double precision      :: etot_save
  double precision      :: si(si_cnt)
  double precision      :: sit(si_cnt), sit_tmp(si_cnt)
  double precision      :: sit2(si_cnt), sit2_tmp(si_cnt)

#ifdef MPI
  double precision              :: new_list_cnt
  double precision, save        :: reduce_buf_in(7)
  double precision, save        :: reduce_buf_out(7)
#endif
  double precision      :: fac(3)
  double precision      :: rmu(3)
  double precision      :: factt
  double precision      :: pconv
  double precision      :: rndf, rndfp, rndfs
  double precision      :: scaltp
  double precision      :: tmassinv
  double precision      :: ocm(3), xcm(3), vcm(3)       ! for COM velocities
  double precision      :: tspan
  double precision      :: velocity2
  double precision      :: vmax
  double precision      :: wfac
  double precision      :: winf
  double precision      :: sys_x, sys_y, sys_z
  double precision      :: sys_range(3, 2)

#ifdef MPI
  integer               :: atm_lst_idx
  integer               :: buf_ctr
#endif
  integer               :: i, j, m
  integer               :: nrx
  integer               :: nstep
  integer               :: steps_per_avg
  integer               :: nvalid
  integer               :: nvalidi

  logical               :: belly
  logical               :: lout
  logical               :: new_list
  logical               :: onstep
  logical               :: reset_velocities
  logical               :: use_vlimit
  logical               :: vlimit_exceeded
  logical               :: write_restrt
  logical               :: doing_langevin
  logical               :: need_pot_enes        ! refers to potential energy...
  logical               :: need_virials

#ifdef MPI
  logical               :: collect_crds
  logical               :: collect_vels
  logical               :: all_crds_valid
  logical               :: all_vels_valid
#endif

  double precision, parameter   :: one_third = 1.d0/3.d0
  double precision, parameter   :: pressure_constant = 6.85695d+4

! xyz components of pressure, ekcmt, virial:

  double precision      :: ekcmt(3)
  double precision      :: press(3)
  double precision      :: virial(3)

! Variables used for Langevin dynamics:

  double precision      :: c_ave
  double precision      :: gammai

! Runmd operates in kcal/mol units for energy, amu for masses,
! and angstoms for distances.  To convert the input time parameters 
! from picoseconds to internal units, multiply by 20.455 
! (which is 10.0 * sqrt(4.184)).

  call zero_time()

! Initialize some variables:

  ! We init for case where they are not used but may be referenced in mpi.

  eke = 0.d0
  ekph = 0.d0
  ekpbs = 0.d0
  ekcmt(:) = 0.d0
  press(:) = 0.d0       ! for ntp 0, this prevents dumping of junk.
  virial(:) = 0.d0

  belly = ibelly .gt. 0
  lout = .true.
  use_vlimit = vlimit .gt. 1.0d-7
  ekmh = 0.d0

  tmassinv = 1.d0 / tmass       ! total mass inverse for all atoms.

! If ntwprt .ne. 0, only print the solute atoms in the coordinate/vel. archives.

  nrx  = natom * 3
  if (ntwprt .gt. 0) nrx = ntwprt * 3

! Get the degrees of freedom for the solute (rndfp) and solvent (rndfs).
! Then correct the solute value for the minimum degrees of freedom in the
! system (due to removal of translational or rotational motion, etc.).

  call degcnt(ibelly, natom, atm_igroup, natom, gbl_bond, gbl_bond(bonda_idx), &
              ntc, rndfp, rndfs)

  rndfp = rndfp - dble(ndfmin)

  rndf = rndfp + rndfs  ! total degrees of freedom in system

! Correct the totaldegrees of freedom for extra points (lone pairs):

  rndf = rndf - 3.d0 * dble(numextra)

! BUGBUG - NOTE that rndfp, rndfs are uncorrected in an extra points context!

! Runmd operates in kcal/mol units for energy, amu for masses,
! and angstoms for distances.  To convert the input time parameters
! from picoseconds to internal units, multiply by 20.455
! (which is 10.0 * sqrt(4.184)).

  boltz2 = 8.31441d-3 * 0.5d0 / 4.184d0

  pconv = 1.6604345d+04
  pconv = pconv * 4.184d0 ! factor to convert pressure from kcal/mole to bar.

  dtx = dt * 20.455d+00
  dtx_inv = 1.0d0 / dtx
  half_dtx = dtx * 0.5d0

  fac(1) = boltz2 * rndf
  fac(2) = boltz2 * rndfp
  if (rndfp .lt. 0.1d0) fac(2) = 1.d-6
  fac(3) = boltz2 * rndfs
  if (rndfs .lt. 0.1d0) fac(3) = 1.d-6

  factt = rndf / (rndf + ndfmin)

  ekin0  = fac(1) * temp0

! Langevin dynamics setup:

  doing_langevin = (gamma_ln .gt. 0.d0)

  gammai = gamma_ln / 20.455d0
  c_ave    = 1.d0 + gammai * half_dtx

  if (doing_langevin .and. ifbox .eq. 0) &
    call get_position(atm_cnt, crd, sys_x, sys_y, sys_z, sys_range)

  if (ntp .gt. 0) dtcp = comp * 1.0d-06 * dt / taup

  nvalid = 0
  nvalidi = 0
  nstep = 0
  steps_per_avg = ene_avg_sampling

  si(:) = 0.d0
  sit(:) = 0.d0
  sit2(:) = 0.d0

! The _tmp variables below are only used if ntave .gt. 0.  They are used for
! scratch space in the ntave calcs, and to hold the last sit* values between
! calls to the ntave code.

  sit_tmp(:) = 0.d0     
  sit2_tmp(:) = 0.d0

#ifdef MPI
! The following flags keep track of coordinate/velocity update status when
! using mpi

  all_crds_valid = .true.
  all_vels_valid = .true.
#endif

! Clean up the velocity if belly run:

  if (belly) call bellyf(atm_cnt, atm_igroup, vel)

! Make a first dynamics step:

  irespa = 1    ! PME respa step counter.

  ! We must build a new pairlist the first time we run force:

  new_list = .true.
  need_pot_enes = .true.        ! maybe not, but most simple for step 0
  need_virials = (ntp .gt. 0)

! irest = 0:  General startup if not continuing a previous run:


  if (irest .eq. 0) then

    irespa = 0

    ! Calculate the force.  This also does ekcmt if a regular pme run:

    call update_time(runmd_time)

    if (using_pme_potential) then

#ifdef MPI
      call pme_force(atm_cnt, crd, frc, gbl_img_atm_map, gbl_atm_img_map, &
                     my_atm_lst, new_list, need_pot_enes, need_virials, &
                     pme_pot_ene, virial, ekcmt, pme_err_est)
#else
      call pme_force(atm_cnt, crd, frc, gbl_img_atm_map, gbl_atm_img_map, &
                     new_list, need_pot_enes, need_virials, &
                     pme_pot_ene, virial, ekcmt, pme_err_est)
#endif /* MPI */

      ! Store energy terms in state info array for printout.
  
      si(si_pot_ene) = pme_pot_ene%total
      si(si_vdw_ene) = pme_pot_ene%vdw_tot
      si(si_elect_ene) = pme_pot_ene%elec_tot
      si(si_hbond_ene) = pme_pot_ene%hbond
      si(si_bond_ene) = pme_pot_ene%bond
      si(si_angle_ene) = pme_pot_ene%angle
      si(si_dihedral_ene) = pme_pot_ene%dihedral
      si(si_vdw_14_ene) = pme_pot_ene%vdw_14
      si(si_elect_14_ene) = pme_pot_ene%elec_14
      si(si_restraint_ene) = pme_pot_ene%restraint
      si(si_pme_err_est) = pme_err_est

      ! Total virial will be printed in mden regardless of ntp value.
      si(si_tot_virial) = virial(1) + virial(2) + virial(3)

    else if (using_gb_potential) then

      call gb_force(atm_cnt, crd, frc, gb_pot_ene, irespa)

      si(si_pot_ene) = gb_pot_ene%total
      si(si_vdw_ene) = gb_pot_ene%vdw_tot
      si(si_elect_ene) = gb_pot_ene%elec_tot
      si(si_hbond_ene) = gb_pot_ene%gb                ! temporary hack
      si(si_bond_ene) = gb_pot_ene%bond
      si(si_angle_ene) = gb_pot_ene%angle
      si(si_dihedral_ene) = gb_pot_ene%dihedral
      si(si_vdw_14_ene) = gb_pot_ene%vdw_14
      si(si_elect_14_ene) = gb_pot_ene%elec_14
      si(si_restraint_ene) = gb_pot_ene%restraint
      si(si_pme_err_est) = 0.d0

    end if

    ! This force call does not count as a "step". call nmrdcp to decrement
    ! local NMR step counter:

    call nmrdcp

    irespa = 1

    ! The coordinates will not be changed between here and the next
    ! run of force, so we can just set new_list to .false.

    new_list = .false.

    ! Reset quantities depending on temp0 and tautp (which may have been 
    ! changed by modwt during force call):

    ekin0  = fac(1) * temp0

    if (ntp .gt. 0) then

      si(si_volume) = uc_volume
      si(si_density) = tmass / (0.602204d0 * si(si_volume))
      si(si_tot_ekcmt) = 0.d0
      si(si_tot_press) = 0.d0

      do m = 1, 3
        ekcmt(m) = ekcmt(m) * 0.5d0
        si(si_tot_ekcmt) = si(si_tot_ekcmt) + ekcmt(m)
        press(m) = (pconv + pconv) * (ekcmt(m) - virial(m)) / si(si_volume)
        si(si_tot_press) = si(si_tot_press) + press(m)
      end do

      si(si_tot_press) = si(si_tot_press) / 3.d0

    end if

    eke = 0.d0

#ifdef MPI
    do atm_lst_idx = 1, my_atm_cnt
      j = my_atm_lst(atm_lst_idx)
#else
    do j = 1, atm_cnt
#endif
      eke = eke + mass(j) * &
             (vel(1,j) * vel(1,j) + vel(2,j) * vel(2,j) + vel(3,j) * vel(3,j))
    end do

    eke = 0.5d0 * eke

    ! New addition - is this right?

#ifdef MPI
    call update_time(runmd_time)
! Sum up the partial kinetic energies. The all_reduce is expensive, but only
! occurs once in an irest .eq. 0 run.

    reduce_buf_in(1) = eke
    call mpi_allreduce(reduce_buf_in, reduce_buf_out, 1, mpi_double_precision, &
                       mpi_sum, mpi_comm_world, err_code_mpi)

    eke = reduce_buf_out(1)
    call update_time(fcve_dist_time)
#endif

    si(si_solute_kin_ene) = eke
    si(si_kin_ene) = si(si_solute_kin_ene)
    si(si_tot_ene) = si(si_kin_ene) + si(si_pot_ene)

    eke = 0.d0

#ifdef MPI
    do atm_lst_idx = 1, my_atm_cnt
      j = my_atm_lst(atm_lst_idx)
#else
    do j = 1, atm_cnt
#endif
      winf = atm_mass_inv(j) * half_dtx
      do i = 1, 3
        vel(i, j) = vel(i, j) - frc(i, j) * winf
        if (use_vlimit) vel(i, j) = sign(min(abs(vel(i, j)), vlimit), vel(i, j))
      end do
    end do

    if (ntt .eq. 1) then
      ekmh = max(si(si_solute_kin_ene), fac(1) * 10.d0)
    end if

  end if

! irest = 1:  Continuation of a previous trajectory:
!
! Note: if the last printed energy from the previous trajectory was
!       at time "t", then the restrt file has velocities at time
!       t + 0.5dt, and coordinates at time t + dt

  ekmh = 0.0d0

#ifdef MPI
  do atm_lst_idx = 1, my_atm_cnt
    i = my_atm_lst(atm_lst_idx)
#else
  do i = 1, atm_cnt
#endif
    ekmh = ekmh + mass(i) * &
           (vel(1,i) * vel(1,i) + vel(2,i) * vel(2,i) + vel(3,i) * vel(3,i))
  end do

#ifdef MPI
! Sum up ekmh. The all_reduce is expensive, but only occurs once per run.

  call update_time(runmd_time)
  reduce_buf_in(1) = ekmh

  call mpi_allreduce(reduce_buf_in, reduce_buf_out, 1, mpi_double_precision, &
                     mpi_sum, mpi_comm_world, err_code_mpi)

  ekmh = reduce_buf_out(1)
  call update_time(fcve_dist_time)
#endif

  ekmh = ekmh * 0.5d0

  last_vel(:,:) = vel(:,:)
  
  if (irest .eq. 0) then

    ! Print the initial energies and temperatures 

    if (nstep .le. 0 .and. master) then
      call prntmd(nstep, t, si, fac, 7, .false.)
      if (nmropt .ne. 0) call nmrptx(6)
    end if
    if (nstlim .eq. 0) return

  end if

  !=======================================================================
  ! MAIN LOOP FOR PERFORMING THE DYNAMICS STEP:
  ! At this point, the coordinates are a half-step "ahead" of the velocities;
  ! the variable EKMH holds the kinetic energy at these "-1/2" velocities,
  ! which are stored in the array last_vel.
  !=======================================================================

  do 

! Calculate the force. This also does ekcmt if a regular pme run:

! Full energies are only calculated every nrespa steps, and may actually be
! calculated less frequently, depending on needs.

    onstep = mod(irespa, nrespa) .eq. 0

    if (vrand .ne. 0 .and. ntt .eq. 2) then
      need_pot_enes = .true.
    else if (ntwe .gt. 0) then
      need_pot_enes = (mod(nstep+1, steps_per_avg) .eq. 0) .or. &
                      (mod(nstep+1, ntwe) .eq. 0 .and. onstep)
    else
      need_pot_enes = mod(nstep+1, steps_per_avg) .eq. 0
    end if

    call update_time(runmd_time)

    if (using_pme_potential) then

#ifdef MPI
      call pme_force(atm_cnt, crd, frc, gbl_img_atm_map, gbl_atm_img_map, &
                     my_atm_lst, new_list, need_pot_enes, need_virials, &
                     pme_pot_ene, virial, ekcmt, pme_err_est)
#else
      call pme_force(atm_cnt, crd, frc, gbl_img_atm_map, gbl_atm_img_map, &
                     new_list, need_pot_enes, need_virials, &
                     pme_pot_ene, virial, ekcmt, pme_err_est)
#endif /* MPI */

      ! Store energy terms in state info array for printout.
  
      si(si_pot_ene) = pme_pot_ene%total
      si(si_vdw_ene) = pme_pot_ene%vdw_tot
      si(si_elect_ene) = pme_pot_ene%elec_tot
      si(si_hbond_ene) = pme_pot_ene%hbond
      si(si_bond_ene) = pme_pot_ene%bond
      si(si_angle_ene) = pme_pot_ene%angle
      si(si_dihedral_ene) = pme_pot_ene%dihedral
      si(si_vdw_14_ene) = pme_pot_ene%vdw_14
      si(si_elect_14_ene) = pme_pot_ene%elec_14
      si(si_restraint_ene) = pme_pot_ene%restraint
      si(si_pme_err_est) = pme_err_est

      ! Total virial will be printed in mden regardless of ntp value.
      si(si_tot_virial) = virial(1) + virial(2) + virial(3)

    else if (using_gb_potential) then

      call gb_force(atm_cnt, crd, frc, gb_pot_ene, irespa)

      si(si_pot_ene) = gb_pot_ene%total
      si(si_vdw_ene) = gb_pot_ene%vdw_tot
      si(si_elect_ene) = gb_pot_ene%elec_tot
      si(si_hbond_ene) = gb_pot_ene%gb                ! temporary hack
      si(si_bond_ene) = gb_pot_ene%bond
      si(si_angle_ene) = gb_pot_ene%angle
      si(si_dihedral_ene) = gb_pot_ene%dihedral
      si(si_vdw_14_ene) = gb_pot_ene%vdw_14
      si(si_elect_14_ene) = gb_pot_ene%elec_14
      si(si_restraint_ene) = gb_pot_ene%restraint
      si(si_pme_err_est) = 0.d0

    end if

! Reset quantities depending on temp0 and tautp (which may have been 
! changed by modwt during force call):

    ekin0 = fac(1) * temp0

! Pressure coupling:

    if (ntp .gt. 0) then

      si(si_volume) = uc_volume
      si(si_density) = tmass / (0.602204d0 * si(si_volume))
      si(si_tot_ekcmt) = 0.d0
      si(si_tot_press) = 0.d0

      do m = 1, 3
        ekcmt(m) = ekcmt(m) * 0.5d0
        si(si_tot_ekcmt) = si(si_tot_ekcmt) + ekcmt(m)
        press(m) = (pconv + pconv) * (ekcmt(m) - virial(m)) / si(si_volume)
        si(si_tot_press) = si(si_tot_press) + press(m)
      end do

      si(si_tot_press) = si(si_tot_press) / 3.d0
      
    end if

    ! Do randomization of velocities, if needed:
  
    ! Assign new random velocities every Vrand steps, if ntt .eq. 2:

    reset_velocities = .false.

    if (vrand .ne. 0 .and. ntt .eq. 2) then
      if (mod((nstep+1), vrand) .eq.  0) reset_velocities = .true.
    end if

    if (reset_velocities) then
      if (master) then
        write(mdout,'(a,i8)') 'Setting new random velocities at step ', &
                              nstep + 1
      end if
      call vrand_set_velocities(atm_cnt, vel, atm_mass_inv, temp0 * factt)
      if (belly) call bellyf(atm_cnt, atm_igroup, vel)

      ! At this point in the code, the velocities lag the positions
      ! by half a timestep.  If we intend for the velocities to be drawn
      ! from a Maxwell distribution at the timepoint where the positions and
      ! velocities are synchronized, we have to correct these newly
      ! redrawn velocities by backing them up half a step using the
      ! current force.
      ! Note that this fix only works for Newtonian dynamics.

      if (.not. doing_langevin) then
#ifdef MPI
        do atm_lst_idx = 1, my_atm_cnt
          j = my_atm_lst(atm_lst_idx)
#else
        do j = 1, atm_cnt
#endif
          wfac = atm_mass_inv(j) * half_dtx
          vel(:, j) = vel(:, j) - frc(:, j) * wfac
        end do
      end if

    end if  ! (reset_velocities)

! Do the velocity update:
   
    if (.not. doing_langevin) then
      
      ! ---Newtonian dynamics:
      
#ifdef MPI
      do atm_lst_idx = 1, my_atm_cnt
        j = my_atm_lst(atm_lst_idx)
#else
      do j = 1, atm_cnt
#endif
        wfac = atm_mass_inv(j) * dtx
        vel(:, j) = vel(:, j) + frc(:, j) * wfac
      end do
      
    else  !  ntt .eq. 3, we are doing langevin dynamics
      
      ! Simple model for Langevin dynamics, basically taken from
      ! Loncharich, Brooks and Pastor, Biopolymers 32:523-535 (1992),
      ! Eq. 11.  (Note that the first term on the rhs of Eq. 11b
      ! should not be there.)
      
      call langevin_setvel(atm_cnt, vel, frc, mass, atm_mass_inv, &
                           dt, temp0, gamma_ln)

    end if      ! (doing_langevin)

! Consider vlimit:

    if (use_vlimit) then

      ! We here provide code that is most efficient if vlimit is not exceeded.

      vlimit_exceeded = .false.

#ifdef MPI
      do atm_lst_idx = 1, my_atm_cnt
        j = my_atm_lst(atm_lst_idx)
#else
      do j = 1, atm_cnt
#endif
        if (abs(vel(1, j)) .gt. vlimit .or. &
            abs(vel(2, j)) .gt. vlimit .or. &
            abs(vel(3, j)) .gt. vlimit) then
          vlimit_exceeded = .true.
          exit
        end if
      end do

      if (vlimit_exceeded) then
        vmax = 0.d0
#ifdef MPI
        do atm_lst_idx = 1, my_atm_cnt
          j = my_atm_lst(atm_lst_idx)
#else
        do j = 1, atm_cnt
#endif
          do i = 1, 3
            vmax = max(vmax, abs(vel(i, j)))
            vel(i, j) = sign(min(abs(vel(i, j)), vlimit), vel(i, j))
          end do
        end do
        ! Only violations on the master node are actually reported
        ! to avoid both MPI communication and non-master writes.
        write(mdout, '(a,i6,a,f10.4)')  'vlimit exceeded for step ', nstep, &
                                        '; vmax = ', vmax
      end if

    end if

! Update the positions, putting the "old" positions into frc():

#ifdef MPI
    do atm_lst_idx = 1, my_atm_cnt
      i = my_atm_lst(atm_lst_idx)
#else
    do i = 1, atm_cnt
#endif
      frc(:, i) = crd(:, i)
      crd(:, i) = crd(:, i) + vel(:, i) * dtx
    end do

! If shake is being used, update new positions to fix bond lengths:

    if (ntc .ne. 1) then

      call update_time(runmd_time)
      call shake(frc, crd)
      call shake_fastwater(frc, crd)
      call update_time(shake_time)

      ! Must update extra point coordinates after the frame is moved:

      if (numextra .gt. 0 ) then
        if (frameon .ne. 0 .and. gbl_frame_cnt .ne. 0) &
          call local_to_global(crd, ep_frames, ep_lcl_crd, gbl_frame_cnt)
      end if

      ! Re-estimate velocities from differences in positions:

#ifdef MPI
      do atm_lst_idx = 1, my_atm_cnt
        j = my_atm_lst(atm_lst_idx)
#else
      do j = 1, atm_cnt
#endif
        vel(:, j) = (crd(:, j) - frc(:, j)) * dtx_inv
      end do

    else

      ! Must update extra point coordinates after the frame is moved:

      if (numextra .gt. 0 ) then
        if (frameon .ne. 0 .and. gbl_frame_cnt .ne. 0) &
          call local_to_global(crd, ep_frames, ep_lcl_crd, gbl_frame_cnt)
      end if

    end if

    ! Extra points velocities are known to be bogus at this point.  We zero
    ! them when they are output; ignore them otherwise.

    if (ntt .eq. 1 .or. onstep) then

      ! Get the kinetic energy, either for printing or for Berendsen.
      ! The process is completed later under mpi, in order to avoid an
      ! extra all_reduce.

      eke = 0.d0
      ekph = 0.d0
      ekpbs = 0.d0

      if (.not. doing_langevin) then

#ifdef MPI
        do atm_lst_idx = 1, my_atm_cnt
          j = my_atm_lst(atm_lst_idx)
#else
        do j = 1, atm_cnt
#endif
          eke = eke + mass(j) * &
                 ((vel(1,j) + last_vel(1,j))**2 + &
                  (vel(2,j) + last_vel(2,j))**2 + &
                  (vel(3,j) + last_vel(3,j))**2)
          
          ekpbs = ekpbs + mass(j) * &
                 ((vel(1,j) * last_vel(1,j)) + &
                  (vel(2,j) * last_vel(2,j)) + &
                  (vel(3,j) * last_vel(3,j)))
          
          ekph  = ekph  + mass(j) * &
                 ((vel(1,j) * vel(1,j)) + &
                  (vel(2,j) * vel(2,j)) + &
                  (vel(3,j) * vel(3,j)))
          
        end do

        eke = 0.25d0 * eke

      else
#ifdef MPI
        do atm_lst_idx = 1, my_atm_cnt
          j = my_atm_lst(atm_lst_idx)
#else
        do j = 1, atm_cnt
#endif
          eke = eke + mass(j) * 0.25d0 * c_ave *&
                 ((vel(1,j) + last_vel(1,j))**2 + &
                  (vel(2,j) + last_vel(2,j))**2 + &
                  (vel(3,j) + last_vel(3,j))**2)
        end do

      end if ! (doing_langevin)
  
      eke = 0.5d0 * eke
      ekph = ekph * 0.5d0
      ekpbs = ekpbs * 0.5d0

!     write(0,*)'eke,ekph,ekpbs=', eke,ekph,ekpbs       ! DBG

      ! NOTE - These ke terms are not yet summed under MPI!!!

    end if ! (ntt .eq. 1 .or. onstep)

! Scale coordinates if constant pressure run:

    if (ntp .eq. 1) then
      rmu(1) = (1.d0 - dtcp * (pres0 - si(si_tot_press)))**one_third
      rmu(2) = rmu(1)
      rmu(3) = rmu(1)
    else if (ntp .gt. 1) then
      rmu(1) = (1.d0 - dtcp * (pres0 - press(1)))**one_third
      rmu(2) = (1.d0 - dtcp * (pres0 - press(2)))**one_third
      rmu(3) = (1.d0 - dtcp * (pres0 - press(3)))**one_third
    end if

    if (ntp .gt. 0) then

      ! WARNING!!   This is not correct for non-orthogonal boxes if NTP > 1
      ! (i.e. non-isotropic scaling).  Currently general cell updates which
      ! allow cell angles to change are not implemented.  The virial tensor
      ! computed for ewald is the general Nose Klein; however the cell 
      ! response needs a more general treatment.

      call pressure_scale_pbc_data(rmu, vdw_cutoff + skinnb, verbose)
      call set_cit_tbl_dims(pbc_box, vdw_cutoff + skinnb, cut_factor)

#ifdef MPI
      call pressure_scale_crds(crd, mass, gbl_mol_mass_inv, gbl_my_mol_lst, &
                               gbl_mol_com)
#else      
      call pressure_scale_crds(crd, mass, gbl_mol_mass_inv, gbl_mol_com)
#endif /* MPI */

      if (ntr .gt. 0 .and. natc .gt. 0) then 
#ifdef MPI
        call pressure_scale_restraint_crds(atm_xc, mass, gbl_mol_mass_inv, &
                                             gbl_my_mol_lst)
#else      
        call pressure_scale_restraint_crds(atm_xc, mass, gbl_mol_mass_inv)
#endif /* MPI */
      end if

    end if        ! ntp .gt. 0

    if (using_pme_potential) then

      ! Now we can do a skin check to see if we will have to rebuild the
      ! pairlist next time...

#ifdef MPI
      call check_my_atom_movement(crd, gbl_atm_saved_crd, my_atm_lst, &
                                  skinnb, ntp, new_list)
#else
      call check_all_atom_movement(atm_cnt, crd, gbl_atm_saved_crd, skinnb, &
                                   ntp, new_list)
#endif /* MPI */

    end if
    ! Nothing to do for Generalized Born...

#ifdef MPI
    if (new_list) then
      new_list_cnt = 1.d0
    else
      new_list_cnt = 0.d0
    end if

! Sum up the partial kinetic energies and complete the skin check. 

    if (ntt .eq. 1 .or. onstep) then
      reduce_buf_in(1) = eke
      reduce_buf_in(2) = ekph
      reduce_buf_in(3) = ekpbs
      reduce_buf_in(4) = new_list_cnt
      buf_ctr = 4
    else
      reduce_buf_in(1) = new_list_cnt
      buf_ctr = 1
    end if

    call update_time(runmd_time)
#ifdef COMM_TIME_TEST
  call start_test_timer(7, 'allreduce eke,list_cnt', 0)
#endif

    call mpi_allreduce(reduce_buf_in, reduce_buf_out, buf_ctr, &
                       mpi_double_precision, &
                       mpi_sum, mpi_comm_world, err_code_mpi)

#ifdef COMM_TIME_TEST
  call stop_test_timer(7)
#endif

    call update_time(fcve_dist_time)

    if (ntt .eq. 1 .or. onstep) then
      eke = reduce_buf_out(1)
      ekph = reduce_buf_out(2)
      ekpbs = reduce_buf_out(3)
      new_list_cnt = reduce_buf_out(4)
    else
      new_list_cnt = reduce_buf_out(1)
    end if

    ! Determine if any process saw an atom exceed the skin check.  We use the
    ! comparison to 0.5d0 to handle any rounding error issues; we use
    ! double precision for new_list_cnt in order to be able to piggyback the 
    ! reduce.

    new_list = (new_list_cnt .ge. 0.5d0)

    if (using_pme_potential) then

      call check_new_list_limit(new_list)

    end if
    ! Nothing to do for Generalized Born...

! Now distribute the coordinates:

    call update_time(runmd_time)
    if (using_gb_potential) then

      ! Generalized Born has it's own coordinate distribution scheme...

      call gb_mpi_allgathervec(atm_cnt, crd)

      all_crds_valid = .true.

    else if (new_list .or. nmropt .ne. 0) then

      ! This always gets done for a work redistribution cycle:

      call mpi_allgathervec(atm_cnt, crd)

      ! If this is a const pressure run and there are coordinate constraints,
      ! we will need to send the adjusted coordinate constraints around when
      ! we redistribute the work load.  Bummer.

      if (ntp .gt. 0) then
        if (ntr .gt. 0 .and. natc .gt. 0) then 
          if (atm_redist_needed .or. fft_redist_needed) then
            call mpi_allgathervec(atm_cnt, atm_xc)
          end if
        end if
      end if

      all_crds_valid = .true.

    else

      call distribute_crds(natom, crd)

      ! It may be necessary for the master to have a complete copy of the
      ! coordinates for archiving and writing the restart file.

      collect_crds = .false.

      if (ntwx .gt. 0) then
        if (mod(nstep + 1, ntwx) .eq. 0) collect_crds = .true.
      end if

      if (ntwr .ne. 0) then
        if (mod(nstep + 1, ntwr) .eq. 0) collect_crds = .true.
      end if

      if (nstep + 1 .ge. nstlim) collect_crds = .true.

      if (collect_crds) call mpi_gathervec(atm_cnt, crd)

      all_crds_valid = .false.

    end if
    call update_time(fcve_dist_time)

#endif /* MPI */

    if (ntt .eq. 1) then
                            
    !           --- following is from T.E. Cheatham, III and B.R. Brooks,
    !               Theor. Chem. Acc. 99:279, 1998.
                                          
#ifdef MPI
      if (my_atm_cnt .gt. 0) &
        scaltp = sqrt(1.d0 + 2.d0 * (dt / tautp) * (ekin0 - eke)/(ekmh + ekph))
                                                        
      do atm_lst_idx = 1, my_atm_cnt
        j = my_atm_lst(atm_lst_idx)
#else
      scaltp = sqrt(1.d0 + 2.d0 * (dt / tautp) * (ekin0 - eke)/(ekmh + ekph))
      do j = 1, atm_cnt
#endif
        vel(:, j) = vel(:, j) * scaltp
      end do

    end if

    ! Pastor, Brooks, Szabo conserved quantity for harmonic oscillator:
    ! Eqn. 4.7b of Mol. Phys. 65:1409-1419, 1988:

    si(si_solvent_kin_ene) = ekpbs + si(si_pot_ene)
    
    si(si_solute_kin_ene) = eke
    si(si_kin_ene) = si(si_solute_kin_ene)

    if (ntt .eq. 1 .and. onstep) then
      ekmh = max(ekph, fac(1) * 10.d0)
    end if
   
    ! If velocities were reset, the KE is not accurate; fudge it here to keep
    ! the same total energy as on the previous step.  Note that this only
    ! affects printout and averages for Etot and KE -- it has no effect on the
    ! trajectory, or on any averages of potential energy terms.
   
    if (reset_velocities) si(si_kin_ene) = etot_save - si(si_pot_ene)
   
    ! --- total energy is sum of KE + PE:
   
    si(si_tot_ene) = si(si_kin_ene) + si(si_pot_ene)
    etot_save = si(si_tot_ene)

! For periodic PME, Zero COM velocity if requested; used for preventing ewald
! "block of ice flying thru space" phenomenon. We make this correction for pme
! explicitly before collecting velocities to the master...

    if (mod(nstep + 1, nscm) .eq. 0) then
      if (ifbox .ne. 0) then
        if (.not. doing_langevin) then

          vcm(:) = 0.d0

#ifdef MPI
          do atm_lst_idx = 1, my_atm_cnt
            j = my_atm_lst(atm_lst_idx)
#else
          do j = 1, atm_cnt
#endif
            aamass = mass(j)
            vcm(:) = vcm(:) + aamass * vel(:, j)
          end do

#ifdef MPI
          call update_time(runmd_time)
          reduce_buf_in(1:3) = vcm(:)

#ifdef COMM_TIME_TEST
  call start_test_timer(8, 'allreduce vcm', 0)
#endif

          call mpi_allreduce(reduce_buf_in, reduce_buf_out, 3, &
                             mpi_double_precision, &
                             mpi_sum, mpi_comm_world, err_code_mpi)

#ifdef COMM_TIME_TEST
  call stop_test_timer(8)
#endif

          vcm(:) = reduce_buf_out(1:3)
          call update_time(fcve_dist_time)
#endif /* MPI */
          vcm(:) = vcm(:) * tmassinv

          if (master) then
            velocity2 = vcm(1) * vcm(1) + vcm(2) * vcm(2) + vcm(3) * vcm(3)
            write(mdout,'(a,f15.6,f9.2,a)') 'check COM velocity, temp: ', &
                                        sqrt(velocity2), 0.5d0 * tmass * &
                                        velocity2 / fac(1), '(Removed)'
          end if

#ifdef MPI
          do atm_lst_idx = 1, my_atm_cnt
            j = my_atm_lst(atm_lst_idx)
#else
          do j = 1, atm_cnt
#endif
            vel(:, j) = vel(:, j) - vcm(:)
          end do
        end if
      end if
    end if      ! (mod(nstep + 1, nscm) .eq. 0)

#ifdef MPI
! It may be necessary for the master to have a complete copy of the velocities
! for archiving and writing the restart file. In any case, all processors must
! have a complete copy of velocities if atom redistribution is going to happen
! in the next call to force.

    ! For GB runs, always execute the "else"

    call update_time(runmd_time)

    if (using_pme_potential .and. new_list .and. &
        (atm_redist_needed .or. fft_redist_needed)) then
      
      if (numextra .gt. 0 .and. frameon .ne. 0) &
        call zero_extra_pnts_vec(vel, ep_frames, gbl_frame_cnt)

      call mpi_allgathervec(atm_cnt, vel)

      all_vels_valid = .true.

    else

      collect_vels = .false.

      if (ntwv .gt. 0) then
        if (mod(nstep + 1, ntwv) .eq. 0) collect_vels = .true.
      else if (collect_crds .and. ntwv .lt. 0) then
        collect_vels = .true.
      end if

      if (ntwr .ne. 0) then
        if (mod(nstep + 1, ntwr) .eq. 0) collect_vels = .true.
      end if

      if (nstep + 1 .ge. nstlim) collect_vels = .true.

      if (collect_vels) then

        if (numextra .gt. 0 .and. frameon .ne. 0) &
          call zero_extra_pnts_vec(vel, ep_frames, gbl_frame_cnt)

        if (using_pme_potential) then
          call mpi_gathervec(atm_cnt, vel)
        else if (using_gb_potential) then
          call gb_mpi_gathervec(atm_cnt, vel)
        end if

      end if

      all_vels_valid = .false.

    end if

    call update_time(fcve_dist_time)
#endif /* MPI */

! Zero COM velocity if requested; here we are doing the adjustment for a
! nonperiodic system, not pme.

    if (mod(nstep + 1, nscm) .eq. 0) then
      if (ifbox .eq. 0) then
#ifdef MPI
        ! WARNING - currently only GB code has ifbox .eq. 0, and currently
        ! all coordinates are always valid for GB.  We include the conditional
        ! below more-or-less as maintenance insurance...
        if (.not. all_crds_valid) then ! Currently always false...
          call update_time(runmd_time)
          call gb_mpi_allgathervec(atm_cnt, crd)
          all_crds_valid = .true.
          call update_time(fcve_dist_time)
        end if
#endif /* MPI */
        if (.not. doing_langevin) then

#ifdef MPI
          ! WARNING - currently GB code never updates all velocities unless
          !           forced by this scenario...
          if (.not. all_vels_valid) then ! Currently always true...

            ! The following conditional is currently never .true., but we
            ! include it if there is ever extra points support under GB, or
            ! if other types of nonperiodic simulations with extra points
            ! are eventually supported...

            if (numextra .gt. 0 .and. frameon .ne. 0) &
              call zero_extra_pnts_vec(vel, ep_frames, gbl_frame_cnt)

            call update_time(runmd_time)
            call gb_mpi_allgathervec(atm_cnt, vel)
            all_vels_valid = .true.
            call update_time(fcve_dist_time)
          end if
#endif /* MPI */
          ! Nonperiodic simulation.  Remove both translation and rotation.
          ! Back the coords up a half step so they correspond to the velocities,
          ! temporarily storing them in frc(:,:).

          frc(:,:) = crd(:,:) - vel(:,:) * half_dtx

          ! Now compute COM motion and remove it; then recompute (sander
          ! compatibility...).

          ! NOTE - if mass can change, that has to be taken into account for
          !        tmass, tmassinv (say for TI).

          call cenmas(atm_cnt, frc, vel, tmass, tmassinv, mass, xcm, vcm, ocm)
          call stopcm(atm_cnt, frc, vel, xcm, vcm, ocm)
          call cenmas(atm_cnt, frc, vel, tmass, tmassinv, mass, xcm, vcm, ocm)

        else

          ! Get current center of the system:

          call get_position(atm_cnt, crd, vcm(1), vcm(2), vcm(3), sys_range)

          ! Recenter system to the original center:

          call re_position(atm_cnt, 0, crd, crd, vcm(1), vcm(2), vcm(3), &
                           sys_x, sys_y, sys_z, sys_range, 0)

        end if
      end if
    end if      ! (mod(nstep + 1, nscm) .eq. 0)

! Zero out any non-moving velocities if a belly is active:

    if (belly) call bellyf(atm_cnt, atm_igroup, vel)

! Save old, or last velocities...

#ifdef MPI
    if (all_vels_valid) then
      last_vel(:,:) = vel(:,:)
    else
      do atm_lst_idx = 1, my_atm_cnt
        j = my_atm_lst(atm_lst_idx)
        last_vel(:, j) = vel(:, j)
      end do
    end if
#else
    last_vel(:,:) = vel(:,:)
#endif /* MPI */

! Update the step counter and the integration time:

    nstep = nstep + 1
    t = t + dt

#ifdef MPI
  ! Here we update running logging averages as appropriate:

  if (loadbal_verbose .gt. 1) then
  
    log_owned_img_cnt_avg = log_owned_img_cnt_avg + &
                            (log_owned_img_cnt * log_recip_nstep)

    log_used_img_cnt_avg = log_used_img_cnt_avg + &
                            (log_used_img_cnt * log_recip_nstep)

    log_owned_atm_cnt_avg = log_owned_atm_cnt_avg + &
                            (log_owned_atm_cnt * log_recip_nstep)

    log_used_atm_cnt_avg = log_used_atm_cnt_avg + &
                            (log_used_atm_cnt * log_recip_nstep)

    log_used_atm_source_cnt_avg = log_used_atm_source_cnt_avg + &
                            (log_used_atm_source_cnt * log_recip_nstep)

    log_provided_atm_cnt_avg = log_provided_atm_cnt_avg + &
                            (log_provided_atm_cnt * log_recip_nstep)

    log_provided_atm_sink_cnt_avg = log_provided_atm_sink_cnt_avg + &
                            (log_provided_atm_sink_cnt * log_recip_nstep)

  end if
#endif

! nvalid is the number of steps where all energies are calculated; the onstep
! variable indicates that the kinetic energy is available (when using respa).

    if (mod(nstep, steps_per_avg) .eq. 0) then

      nvalid = nvalid + 1

      if (master) then
        do m = 1, si_cnt
          sit(m) = sit(m) + si(m)
          sit2(m) = sit2(m) + si(m) * si(m)
        end do
      end if

    end if

    lout = mod(nstep, ntpr) .eq. 0 .and. onstep

    ! added for rbornstat
    if (mod(irespa,nrespai) == 0 .or. irespa < 2) nvalidi = nvalidi + 1

    irespa = irespa + 1

! Output from this step if required:

! Only the master needs to do the output:

    if (master) then

! Restrt:

      write_restrt = .false.

      if (nstep .eq. nstlim) then
        write_restrt = .true.
      else if (ntwr .ne. 0) then
        if (mod(nstep, ntwr) .eq. 0 .and. onstep) write_restrt = .true.
      end if

      if (write_restrt) then

#ifdef MPI
#else
        if (numextra .gt. 0 .and. frameon .ne. 0) &
          call zero_extra_pnts_vec(vel, ep_frames, gbl_frame_cnt)
#endif

        if (iwrap .eq. 0) then
          call mdwrit(nstep, atm_cnt, crd, vel, t)
        else
          call wrapped_mdwrit(nstep, atm_cnt, crd, vel, t)
        end if
      end if

! Coordinate archive:

      if (ntwx .gt. 0) then

        if (mod(nstep, ntwx) .eq. 0) then

          if (iwrap .eq. 0) then
            call corpac(nrx, crd, 1, mdcrd)
          else
            call wrapped_corpac(atm_cnt, nrx, crd, 1, mdcrd)
          end if

          if (ntb .gt. 0) then
            if (ioutfm .eq. 1) then
              call write_binary_cell_dat
            else
              call corpac(3, pbc_box, 1, mdcrd)
            end if 
          end if

          ! For ioutfm == 1, coordinate archive may also contain vels.

          if (ioutfm .eq. 1 .and. ntwv .lt. 0) then
#ifdef MPI
#else
            ! The zeroing of extra points may be redundant here.
            ! Not a big deal overall though.
            if (numextra .gt. 0 .and. frameon .ne. 0) &
              call zero_extra_pnts_vec(vel, ep_frames, gbl_frame_cnt)
#endif /* MPI */

            call corpac(nrx, vel, 1, mdvel)

          end if

        end if

      end if

! Velocity archive:

      if (ntwv .gt. 0) then

        if (mod(nstep, ntwv) .eq. 0) then
#ifdef MPI
#else
          ! The zeroing of extra points may be redundant here.
          ! Not a big deal overall though.
          if (numextra .gt. 0 .and. frameon .ne. 0) &
            call zero_extra_pnts_vec(vel, ep_frames, gbl_frame_cnt)
#endif /* MPI */

          call corpac(nrx, vel, 1, mdvel)

        end if

      end if

! Energy archive:

      if (ntwe .gt. 0) then
        if (mod(nstep, ntwe) .eq. 0 .and. onstep) then
          call mdeng(nstep, t, si, fac, press, virial, ekcmt)
        end if
      end if

      if (lout) then
        call prntmd(nstep, t, si, fac, 7, .false.)
        if (nmropt .ne. 0) call nmrptx(6)
      end if

! Flush binary netCDF file(s) if necessary...

      if (ioutfm .eq. 1) then

        if (ntwx .gt. 0) then
          if (mod(nstep, ntwx) .eq. 0) then
            call end_binary_frame(mdcrd)
          end if
        end if

        if (ntwv .gt. 0) then

          if (mod(nstep, ntwv) .eq. 0) then
#ifdef MPI
#else
            ! The zeroing of extra points may be redundant here.
            ! Not a big deal overall though.
            if (numextra .gt. 0 .and. frameon .ne. 0) &
              call zero_extra_pnts_vec(vel, ep_frames, gbl_frame_cnt)
#endif /* MPI */

            call end_binary_frame(mdvel)

          end if
        end if

      end if

! Output running average:

      if (ntave .gt. 0) then
        if (mod(nstep, ntave) .eq. 0 .and. onstep) then

          write(mdout, 542)

          tspan = ntave / nrespa

          ! Coming into this loop, the _tmp variables hold the values of
          ! sit, sit2 when this routine was last called (or 0.d0).  The _tmp
          ! vars are used here as scatch space and then updated with the
          ! current sit, sit2.

          do m = 1, si_cnt
            sit_tmp(m) = (sit(m) - sit_tmp(m)) / tspan
            sit2_tmp(m) = (sit2(m) - sit2_tmp(m)) / tspan - &
                          sit_tmp(m) * sit_tmp(m)
            if (sit2_tmp(m) .lt. 0.d0) sit2_tmp(m) = 0.d0
            sit2_tmp(m) = sqrt(sit2_tmp(m))
          end do

          write(mdout, 540) ntave / nrespa
          call prntmd(nstep, t, sit_tmp, fac, 0, .false.)
          write(mdout, 550)
          call prntmd(nstep, t, sit2_tmp, fac, 0, .true.)
          write(mdout, 542)

          sit_tmp(:) = sit(:)
          sit2_tmp(:) = sit2(:)

        end if
      end if

    end if        ! end master's output

    if (nstep .ge. nstlim) exit
    call update_time(runmd_time)

  end do ! Major cycle back to new step unless we have reached our limit:

! Print averages:

  if (master) then

    tspan = nvalid

    if (nvalid .gt. 0) then

      do m = 1, si_cnt
        sit(m) = sit(m)/tspan
        sit2(m) = sit2(m)/tspan - sit(m) * sit(m)
        if (sit2(m) .lt. 0.d0) sit2(m) = 0.d0
        sit2(m) =  sqrt(sit2(m))
      end do

      write(mdout, 540) nvalid
      call prntmd(nstep, t, sit, fac, 0, .false.)
      if (nmropt .ne. 0) call nmrptx(6)
      write(mdout, 550)
      call prntmd(nstep, t, sit2, fac, 0, .true.)

      if (nmropt .ne. 0) then
        write(mdout, '(/,a,/)') ' NMR restraints on final step:'
        call ndvptx(crd, frc, 6)
      end if

      ! Print Born radii statistics:

      if (using_gb_potential) then
        if (rbornstat .eq. 1) then
          ! Born radii stats collected every nrespai step not nrespa step
          tspan = nvalidi
          write(mdout, 580) nstep
          write(mdout, 590)
          do m = 1, atm_cnt
            gbl_rbave(m) = gbl_rbave(m) / tspan
            gbl_rbfluct(m) = gbl_rbfluct(m) / tspan - &
                             gbl_rbave(m) * gbl_rbave(m)
            gbl_rbfluct(m) = sqrt(gbl_rbfluct(m))
            write(mdout, 600) m, gbl_rbmax(m), gbl_rbmin(m), &
                              gbl_rbave(m), gbl_rbfluct(m)
          end do
        end if
      end if

    end if ! (nvalid .gt. 0)

  end if

  call update_time(runmd_time)

  return

  540 format(/5x, ' A V E R A G E S   O V E R ', i7, ' S T E P S', /)
  542 format('|',79('='))
  550 format(/5x, ' R M S  F L U C T U A T I O N S', /)
  580 format('STATISTICS OF EFFECTIVE BORN RADII OVER ',i7,' STEPS')
  590 format('ATOMNUM     MAX RAD     MIN RAD     AVE RAD     FLUCT')
  600 format(i4, 2x, 4f12.4)

end subroutine runmd


!*******************************************************************************
!
! Subroutine:  cenmas
!
! Description:  Calculate the translational and rotational kinetic energies
!               and velocities.
!              
!*******************************************************************************

subroutine cenmas(atm_cnt, x, v, tmass, tmassinv, amass, xcm, vcm, ocm)

  use file_io_dat_mod
  use gbl_constants_mod
  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: x(*)
  double precision      :: v(*)
  double precision      :: tmass
  double precision      :: tmassinv
  double precision      :: amass(*)
  double precision      :: xcm(*)
  double precision      :: vcm(*)
  double precision      :: ocm(*)

! Local variables:
   
  double precision      :: aamass
  double precision      :: comvel
  double precision      :: det
  integer               :: i, j, m, n
  integer               :: i0
  integer               :: lh(3)
  integer               :: mh(3)
  double precision      :: acm(3)
  double precision      :: tcm(3, 3)
  double precision      :: xx, xy, xz, yy, yz, zz
  double precision      :: x1, x2, x3
  double precision      :: ekcm, ekrot

  double precision, parameter   :: crit = 1.d-06
   
  i0 = 3 * atm_cnt
   
  ! calculate the center of mass coordinates:
   
  xcm(1) = 0.d0
  xcm(2) = 0.d0
  xcm(3) = 0.d0
   
  i = 0
  do j = 1, atm_cnt
    aamass = amass(j)
    xcm(1) = xcm(1) + x(i+1) * aamass
    xcm(2) = xcm(2) + x(i+2) * aamass
    xcm(3) = xcm(3) + x(i+3) * aamass
    i = i + 3
  end do

  xcm(1) = xcm(1) * tmassinv
  xcm(2) = xcm(2) * tmassinv
  xcm(3) = xcm(3) * tmassinv
   
  ! calculate velocity and translational kinetic energy of the center of mass:
   
  ekcm = 0.d0
  vcm(1) = 0.0d0
  vcm(2) = 0.0d0
  vcm(3) = 0.0d0
   
  i = 0
  do j = 1, atm_cnt
    aamass = amass(j)
    vcm(1) = vcm(1) + v(i+1) * aamass
    vcm(2) = vcm(2) + v(i+2) * aamass
    vcm(3) = vcm(3) + v(i+3) * aamass
    i = i + 3
  end do
  
  do i = 1, 3
    vcm(i) = vcm(i) * tmassinv
    ekcm = ekcm + vcm(i) * vcm(i)
  end do
  ekcm = ekcm * tmass * 0.5d0
  comvel = sqrt(vcm(1) * vcm(1) + vcm(2) * vcm(2) + vcm(3) * vcm(3))
   
  ! calculate the angular momentum about the center of mass:
   
  acm(1) = 0.0d0
  acm(2) = 0.0d0
  acm(3) = 0.0d0

  i = 0
  do j = 1, atm_cnt
    aamass = amass(j)
    acm(1) = acm(1) + (x(i+2) * v(i+3)-x(i+3) * v(i+2)) * aamass
    acm(2) = acm(2) + (x(i+3) * v(i+1)-x(i+1) * v(i+3)) * aamass
    acm(3) = acm(3) + (x(i+1) * v(i+2)-x(i+2) * v(i+1)) * aamass
    i = i + 3
 end do

  acm(1) = acm(1) - (xcm(2) * vcm(3)-xcm(3) * vcm(2)) * tmass
  acm(2) = acm(2) - (xcm(3) * vcm(1)-xcm(1) * vcm(3)) * tmass
  acm(3) = acm(3) - (xcm(1) * vcm(2)-xcm(2) * vcm(1)) * tmass

  ! calculate the inertia tensor:
   
  xx = 0.d0
  xy = 0.d0
  xz = 0.d0
  yy = 0.d0
  yz = 0.d0
  zz = 0.d0

  i = 0
  do j = 1, atm_cnt
    x1 = x(i+1) - xcm(1)
    x2 = x(i+2) - xcm(2)
    x3 = x(i+3) - xcm(3)
    aamass = amass(j)
    xx = xx + x1 * x1 * aamass
    xy = xy + x1 * x2 * aamass
    xz = xz + x1 * x3 * aamass
    yy = yy + x2 * x2 * aamass
    yz = yz + x2 * x3 * aamass
    zz = zz + x3 * x3 * aamass
    i = i + 3
  end do

  tcm(1,1) = yy + zz
  tcm(2,1) = -xy
  tcm(3,1) = -xz
  tcm(1,2) = -xy
  tcm(2,2) = xx + zz
  tcm(3,2) = -yz
  tcm(1,3) = -xz
  tcm(2,3) = -yz
  tcm(3,3) = xx + yy
   
  ! invert the inertia tensor:
   
  call matinv(tcm, 3, det, lh, mh)

  if (abs(det) .le. crit) then
    write(mdout, '(a,a)') error_hdr, 'Zero determinant calculated in cenmas()'
    call mexit(6, 1)
  end if
   
  ! calculate the angular velocity about the center of mass and the rotational
  ! kinetic energy:
   
  ekrot = 0.d0
  do m = 1, 3
    ocm(m) = 0.d0
    do n = 1, 3
      ocm(m) = ocm(m) + tcm(m,n) * acm(n)
    end do
    ekrot = ekrot + ocm(m) * acm(m)
  end do
  ekrot = ekrot * 0.5d0

  if (master) then
    write(mdout, '(/3x, a, f11.4, 3x, a, f11.4, 3x, a, f12.6)') 'KE Trans =', &
          ekcm, 'KE Rot =', ekrot, 'C.O.M. Vel =', comvel
  end if

  return

end subroutine cenmas 

!*******************************************************************************
!
! Subroutine:   stopcm
!
! Description:  Remove Center of Mass transrotational motion.
!              
!*******************************************************************************

subroutine stopcm(nr, x, v, xcm, vcm, ocm)
   
  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: nr
  double precision      :: x(*)
  double precision      :: v(*)
  double precision      :: xcm(*)
  double precision      :: vcm(*)
  double precision      :: ocm(*)

! Local variables:
   
  integer               :: i, j, m
  double precision      :: x1, x2, x3

  ! stop the center of mass translation:
   
  i = 0
  do j = 1, nr
    do m = 1, 3
      i = i + 1
      v(i) = v(i) - vcm(m)
    end do
  end do
   
  ! stop the rotation about the center of mass:
   
  i = 0
  do j = 1, nr
      x1 = x(i+1) - xcm(1)
      x2 = x(i+2) - xcm(2)
      x3 = x(i+3) - xcm(3)
      v(i+1) = v(i+1) - ocm(2) * x3 + ocm(3) * x2
      v(i+2) = v(i+2) - ocm(3) * x1 + ocm(1) * x3
      v(i+3) = v(i+3) - ocm(1) * x2 + ocm(2) * x1
      i = i + 3
  end do

  if (master) write(mdout, 9008)
9008 format(/3x, 'Translational and rotational motion removed')

  return

end subroutine stopcm 

!*******************************************************************************
!
! Subroutine:   get_position
!
! Description:  Find the center of a set of atoms, and the extent.
!              
! Input:  n  number of atoms
!         x  coordinates
! Output: xc, yc, zc  coordinates geometric center
!         e   xmin, ymin, zmin, xmax, ymax, zmax
!
!*******************************************************************************

subroutine get_position(n, x, xc, yc, zc, e)

  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: n
  double precision      :: x(3,n)
  double precision      :: xc, yc, zc
  double precision      :: e(3,2)

! Local variables:

  integer               :: i, j

  xc = 0.d0
  yc = 0.d0
  zc = 0.d0

  do i = 1, 3
    e(i,1) = x(i,1)
    e(i,2) = x(i,1)
  end do

  do i = 1, n
    do j = 1, 3
      e(j,1) = min(e(j,1), x(j,i))
      e(j,2) = max(e(j,2), x(j,i))
    end do
  end do

  xc = (e(1,2) + e(1,1)) * 0.5d0
  yc = (e(2,2) + e(2,1)) * 0.5d0
  zc = (e(3,2) + e(3,1)) * 0.5d0

  return

end subroutine get_position 

!*******************************************************************************
!
! Subroutine:   re_position
!
! Description:  move the center of a set of atoms
!              
! Input:  n  number of atoms
!         x  coordinates
!         xr reference coordinates
!         xc, yc, zc current center
!         x0, y0, z0 new center
!   
! Output: x, xr  moved coordinates
!   
!*******************************************************************************

subroutine re_position(n, ntr, x, xr, xc, yc, zc, x0, y0, z0, vec, verbose)
   
  use file_io_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: n
  integer               :: ntr
  double precision      :: x(3,n)
  double precision      :: xr(3,n)
  double precision      :: xc, yc, zc
  double precision      :: x0, y0, z0
  double precision      :: vec(3)
  integer               :: verbose

! Local variables:

  integer               :: i, j
  double precision      :: xd, yd, zd

  xd = x0 - xc
  yd = y0 - yc
  zd = z0 - zc

  if (master) &
    write(mdout, '(a, 3f10.6)')"| RE_POSITION Moving by ", xd, yd, zd

  vec(1) = vec(1) + xd
  vec(2) = vec(2) + yd
  vec(3) = vec(3) + zd

  do i = 1, n
    x(1,i) = xd + x(1,i)
    x(2,i) = yd + x(2,i)
    x(3,i) = zd + x(3,i)
  end do

  if (ntr .gt. 0) then
    do i = 1, n
      xr(1,i) = xd + xr(1,i)
      xr(2,i) = yd + xr(2,i)
      xr(3,i) = zd + xr(3,i)
    end do
  end if

  if (verbose .gt. 0 .and. master) then
    write(mdout, *)"*********** Coords moved ****************"
    write(mdout, *) "delta x,y,z ", xd, yd, zd
  end if

  return

end subroutine re_position 

!*******************************************************************************
!
! Subroutine:  matinv
!
! Description: Standard ibm matrix inversion routine.
!
! Arguments:
!   a:  square matrix of dimension nxn
!   d:  resultant determinant
!   l:  work vector of length n
!   m:  work vector of length n
!              
!*******************************************************************************

subroutine matinv(a, n, d, l, m)

  implicit none
!  implicit double precision    :: (a-h,o-z)

! Formal arguments:

  double precision      :: a(*)
  integer               :: n
  double precision      :: d
  integer               :: l(*)
  integer               :: m(*)

! Local variables:

  double precision      :: biga
  double precision      :: hold
  integer               :: i, ij, ik, iz
  integer               :: j, ji, jk, jp, jq, jr
  integer               :: k, ki, kj, kk
  integer               :: nk
   
    ! Search for largest element:
   
    d = 1.0d0
    nk = -n

    do 80 k = 1, n
      nk = nk + n
      l(k) = k
      m(k) = k
      kk = nk + k
      biga = a(kk)
      do 20 j = k, n
        iz = n * (j - 1)
        do 20 i = k, n
          ij = iz + i
          if (abs(biga) - abs(a(ij))) 15, 20, 20
 15         biga = a(ij)
            l(k) = i
            m(k) = j
 20         continue
      
    ! Interchange rows:
      
      j = l(k)
      if (j - k) 35, 35, 25
 25     ki = k - n
        do 30 i = 1, n
          ki = ki + n
          hold = -a(ki)
          ji = ki - k + j
          a(ki) = a(ji)
 30       a(ji) = hold
      
    ! Interchange columns:
      
 35     i = m(k)

      if (i - k) 45, 45, 38
 38     jp = n * (i - 1)
        do 40 j = 1, n
          jk = nk + j
          ji = jp + j
          hold = -a(jk)
          a(jk) = a(ji)
 40       a(ji) = hold
      
    ! Divide column by minus pivot:
      
 45     if (biga) 48, 46, 48
 46       d = 0.0d0
          goto 150
 48       do 55 i = 1, n
            if (i - k) 50, 55, 50
 50           ik = nk + i
              a(ik) = a(ik)/(-biga)
 55           continue
      
     ! Reduce matrix:
      
      do 65 i = 1, n
        ik = nk + i
        hold = a(ik)
        ij = i - n
        do 65 j = 1, n
          ij = ij + n
          if (i - k) 60, 65, 60
 60         if (j - k) 62, 65, 62
 62           kj = ij - i + k
              a(ij) = hold * a(kj) + a(ij)
 65       continue
      
    ! Divide row by pivot:
      
      kj = k - n
      do 75 j = 1, n
        kj = kj + n
        if (j - k) 70, 75, 70
 70       a(kj) = a(kj)/biga
 75     continue
      
    ! Product of pivots:
      
      d = d * biga
      
    ! Replace pivot by reciprocal:
      
      a(kk) = 1.0d0/biga

 80 continue
   
    ! Final row and column interchange:
   
    k = n
100 k = (k - 1)
    if (k) 150, 150, 105
105   i = l(k)
      if (i - k) 120, 120, 108
108     jq = n * (k - 1)
        jr = n * (i - 1)
        do 110 j = 1, n
          jk = jq + j
          hold = a(jk)
          ji = jr + j
          a(jk) = -a(ji)
110       a(ji) = hold
120     j = m(k)
      if (j - k) 100, 100, 125
125     ki = k - n
        do 130 i = 1, n
          ki = ki + n
          hold = a(ki)
          ji = ki - k + j
          a(ki) = -a(ji)
130       a(ji) = hold
    goto 100

150 return

end subroutine matinv 

end module runmd_mod
