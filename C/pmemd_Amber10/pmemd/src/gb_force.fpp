#include "copyright.i"

!*******************************************************************************
!
! Module: gb_force_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module gb_force_mod

  use gbl_datatypes_mod

  implicit none

  ! Potential energies, with breakdown, from GB.  This is intended to be the
  ! external interface to potential energy data produced by this module, in
  ! particular by subroutine gb_force().

  type gb_pot_ene_rec
    sequence
    double precision    :: total
    double precision    :: vdw_tot
    double precision    :: elec_tot
    double precision    :: gb
    double precision    :: bond
    double precision    :: angle
    double precision    :: dihedral
    double precision    :: vdw_14
    double precision    :: elec_14
    double precision    :: restraint
  end type gb_pot_ene_rec

  integer, parameter    :: gb_pot_ene_rec_size = 10

  type(gb_pot_ene_rec), parameter      :: null_gb_pot_ene_rec = &
    gb_pot_ene_rec(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)

contains

!*******************************************************************************
!
! Subroutine:  gb_force
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine gb_force(atm_cnt, crd, frc, pot_ene, irespa)

  use angles_mod
  use bonds_mod
  use constraints_mod
  use dihedrals_mod
  use extra_pnts_nb14_mod
  use gb_ene_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use nmr_calls_mod
  use parallel_dat_mod
  use parallel_mod
  use gb_parallel_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use runfiles_mod
  use timers_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt
  double precision              :: crd(3, atm_cnt)
  double precision              :: frc(3, atm_cnt)
  type(gb_pot_ene_rec)          :: pot_ene
  integer, intent(in)           :: irespa

! Local variables:

  double precision              :: enmr(3)
  integer                       :: i, j

! Zero energies that are stack or call parameters:

  pot_ene = null_gb_pot_ene_rec

! Zero internal energies:

  enmr(:) = 0.d0

! Do weight changes, if requested.

  if (nmropt .ne. 0) call nmr_weight(atm_cnt, crd, 6)

! If no force calcs are to be done, clear the frc array and bag out now.

  if (ntf .eq. 8) then
    frc(:,:) = 0.d0
    return
  end if

  call zero_time()
  call zero_gb_time()

! Calculate the non-bonded contributions:

  frc(:,:) = 0.d0

  call gb_ene(crd, frc, atm_gb_radii, atm_gb_fs, atm_qterm, atm_iac, &
              typ_ico, atm_numex, gbl_natex, atm_cnt, belly_atm_cnt, &
              pot_ene%gb, pot_ene%elec_tot, pot_ene%vdw_tot, irespa)

! Calculate the 1-4 vdw and electrostatics contributions:

  call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
                       gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                       pot_ene%elec_14, pot_ene%vdw_14)

  call update_time(nonbond_time)
  
! Calculate the other contributions:

  call gb_bonded_force(crd, frc, pot_ene)

  ! Sum up total potential energy for this task:

  pot_ene%total = pot_ene%vdw_tot + &
                  pot_ene%elec_tot + &
                  pot_ene%gb + &
                  pot_ene%bond + &
                  pot_ene%angle + &
                  pot_ene%dihedral + &
                  pot_ene%vdw_14 + &
                  pot_ene%elec_14 + &
                  pot_ene%restraint

#ifdef MPI
! Distribute forces to atom owners; this is a reduce sum operation:

  call gb_frcs_distrib(atm_cnt, frc, dbl_mpi_recv_buf)

! Distribute energies to all processes.

  call gb_distribute_enes(pot_ene)

#endif /* MPI */

! Calculate the NMR restraint energy contributions, if requested.
 
  if (nmropt .ne. 0) then
    call nmr_calc(crd, frc, enmr, 6)
    pot_ene%restraint = pot_ene%restraint + enmr(1) + enmr(2) + enmr(3)
    pot_ene%total = pot_ene%total + enmr(1) + enmr(2) + enmr(3)
  end if

! If belly is on then set the belly atom forces to zero:

  if (ibelly .gt. 0) call bellyf(atm_cnt, atm_igroup, frc)

  return

end subroutine gb_force

!*******************************************************************************
!
! Subroutine:  gb_bonded_force
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine gb_bonded_force(crd, frc, pot_ene)

  use angles_mod
  use bonds_mod
  use constraints_mod
  use dihedrals_mod
  use dynamics_dat_mod
  use mdin_ctrl_dat_mod
  use nmr_calls_mod
  use timers_mod

  implicit none

! Formal arguments:

  double precision              :: crd(3, *)
  double precision              :: frc(3, *)
  type(gb_pot_ene_rec)          :: pot_ene

! Local variables:

  ! These energy variables are temporaries, for summing. DON'T use otherwise!

  double precision              :: bond_ene
  double precision              :: angle_ene
  double precision              :: dihedral_ene

  double precision              :: molvir(3, 3)         ! dummy argument
  double precision              :: e14_vir(3, 3)        ! dummy argument

! Bond energy contribution:

! The ebdev/eadev stuff currently only is output under nmr_calls for non-mpi
! code, so we basically drop it here under mpi.

#ifdef MPI
#else
  ebdev = 0.d0
#endif
  if (ntf .le. 1) then
    if (cit_nbonh .gt. 0) then
      call get_bond_energy(cit_nbonh, cit_h_bond, crd, frc, bond_ene)
      pot_ene%bond = bond_ene
    end if
  end if

  if (ntf .le. 2) then
    if (cit_nbona .gt. 0) then
      call get_bond_energy(cit_nbona, cit_a_bond, crd, frc, bond_ene)
      pot_ene%bond = pot_ene%bond + bond_ene
    end if
  end if
#ifdef MPI
#else
    if (cit_nbonh + cit_nbona .gt. 0) &
      ebdev = sqrt(ebdev / (cit_nbonh + cit_nbona))
#endif

  call update_time(bond_time)

! Angle energy contribution:

#ifdef MPI
#else
  eadev = 0.d0
#endif
  if (ntf .le. 3) then
    if (cit_ntheth .gt. 0) then
      call get_angle_energy(cit_ntheth, cit_angle, crd, frc, angle_ene)
      pot_ene%angle = angle_ene
    end if
  end if

  if (ntf .le. 4) then
    if (cit_ntheta .gt. 0) then
      call get_angle_energy(cit_ntheta, cit_angle(cit_ntheth+1), &
                            crd, frc, angle_ene)
      pot_ene%angle = pot_ene%angle + angle_ene
    end if
  end if
#ifdef MPI
#else
  if (cit_ntheth + cit_ntheta .gt. 0) &
    eadev = 57.296 * sqrt(eadev / (cit_ntheth + cit_ntheta))
#endif

  call update_time(angle_time)

! Dihedral energy contribution:

  if (ntf .le. 5) then
    if (cit_nphih .gt. 0) then
      call get_dihed_energy(cit_nphih, cit_dihed, crd, frc, dihedral_ene)
      pot_ene%dihedral = dihedral_ene
    end if
  end if

  if (ntf .le. 6) then
    if (cit_nphia .gt. 0) then
      call get_dihed_energy(cit_nphia, cit_dihed(cit_nphih + 1), crd, frc, &
                            dihedral_ene)
      pot_ene%dihedral = pot_ene%dihedral + dihedral_ene
    end if
  end if

  call update_time(dihedral_time)

! Calculate the position constraint energy:

  if (natc .gt. 0) then
     call get_crd_constraint_energy(natc, pot_ene%restraint, atm_jrc, &
                                    crd, frc, atm_xc, atm_weight)
  endif

  return

end subroutine gb_bonded_force

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  gb_distribute_enes
!
! Description: <TBS>
!
!*******************************************************************************

subroutine gb_distribute_enes(pot_ene)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  type(gb_pot_ene_rec) :: pot_ene

! Local variables:

  type(gb_pot_ene_rec), save    :: dat_in, dat_out

  dat_in = pot_ene

  call mpi_allreduce(dat_in%total, dat_out%total, &
                     gb_pot_ene_rec_size, mpi_double_precision, &
                     mpi_sum, mpi_comm_world, err_code_mpi)

  pot_ene = dat_out

  return

end subroutine gb_distribute_enes
#endif


end module gb_force_mod
