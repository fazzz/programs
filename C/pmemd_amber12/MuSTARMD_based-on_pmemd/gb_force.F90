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
    double precision    :: surf
    double precision    :: bond
    double precision    :: angle
    double precision    :: dihedral
    double precision    :: vdw_14
    double precision    :: elec_14
    double precision    :: restraint
    double precision    :: angle_ub !CHARMM Urey-Bradley Energy
    double precision    :: imp      !CHARMM Improper Energy
    double precision    :: cmap     !CHARMM CMAP

  end type gb_pot_ene_rec

  integer, parameter    :: gb_pot_ene_rec_size = 14

  type(gb_pot_ene_rec), parameter      :: null_gb_pot_ene_rec = &
    gb_pot_ene_rec(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                   0.d0,0.d0,0.d0,0.d0)

contains

!*******************************************************************************
!
! Subroutine:  gb_force
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine gb_force(atm_cnt, crd, frc, pot_ene, irespa, need_pot_ene)

  use angles_mod
  use angles_ub_mod
  use bonds_mod
  use constraints_mod
  use dihedrals_mod
  use dihedrals_imp_mod
  use cmap_mod
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
  use charmm_mod, only : charmm_active
  use mdin_debugf_dat_mod, only : do_charmm_dump_gold
  use amd_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt
  double precision              :: crd(3, atm_cnt)
  double precision              :: frc(3, atm_cnt)
  type(gb_pot_ene_rec)          :: pot_ene
  integer, intent(in)           :: irespa
  logical, intent(in)           :: need_pot_ene

! Local variables:

  double precision              :: enmr(3)
  integer                       :: i, j
  double precision              :: temp_tot_dih
  integer                       :: buf_size


#ifdef CUDA
  ! Update any weight modifications and zero-out the potential energy record
  if (nmropt .ne. 0) then
    call nmr_weight(atm_cnt, crd, 6)
    pot_ene = null_gb_pot_ene_rec !Zero's entire structure
    enmr(:) = 0
  end if
  if ((mod(irespa, ene_avg_sampling) .eq. 0) .or. (mod(irespa, ntpr) .eq. 0) &
      .or. need_pot_ene .or. (iamd.gt.0)) then
#ifdef MPI    
    if((iamd.eq.2).or.(iamd.eq.3))then
      totdih = 0.0
      if(num_amd_lag .eq. 0)then
        temp_tot_dih = 0.0
        !Note this double calculates the dihedral energy since we do it
        !here for AMD and then repeat it later - consider merging these
        !at some point for efficiency.
        call gpu_calculate_gb_amd_dihedral_energy(totdih)
        buf_size = 1
        call mpi_allreduce(totdih, temp_tot_dih, &
                       buf_size, mpi_double_precision, &
                       mpi_sum, pmemd_comm, err_code_mpi)
        totdih = temp_tot_dih
      end if
      call gpu_calculate_amd_dihedral_weight(totdih)
    end if
#else    
    if((iamd.eq.2).or.(iamd.eq.3))then
      call gpu_calculate_gb_amd_dihedral_energy_weight()
    endif
#endif    
    call gpu_gb_ene(pot_ene)
#ifdef MPI    
    call gb_distribute_enes(pot_ene)
#endif    
  else
    call gpu_gb_forces()
  end if
  ! Calculate NMR restraints and upload the forces back to the GPU
  ! Note we ultimately need a GPU version of nmr_calc to avoid these
  ! costly transfers on every step.
  ! Do this every step even if using RESPA.
  if (nmropt .ne. 0) then
    call gpu_download_crd(crd)
    call gpu_download_frc(frc)
    call nmr_calc(crd, frc, enmr, 6)
    ! May not need pot_ene%restraint and total if not printing energies, but
    ! doesn't hurt and avoids additional if statement.
    pot_ene%restraint = pot_ene%restraint + enmr(1) + enmr(2) + enmr(3)
    pot_ene%total = pot_ene%total + enmr(1) + enmr(2) + enmr(3)
    call gpu_upload_frc(frc)
  end if
  if(iamd.gt.0)then
!AMD calculate weight and scale forces
    call gpu_calculate_and_apply_amd_weights(pot_ene%total, pot_ene%dihedral, &
                                             num_amd_lag)
  endif
  call update_time(nonbond_time)
  return
#else

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
              pot_ene%gb, pot_ene%elec_tot, pot_ene%vdw_tot,         &
              pot_ene%surf, irespa)

! Calculate the 1-4 vdw and electrostatics contributions:

  if (charmm_active) then
    call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
                         gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                         pot_ene%elec_14, pot_ene%vdw_14)
  else
    call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
                         gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                         pot_ene%elec_14, pot_ene%vdw_14)
  endif 

  call update_time(nonbond_time)
  
! Calculate the other contributions:

  call gb_bonded_force(crd, frc, pot_ene)

  ! Sum up total potential energy for this task:

  pot_ene%total = pot_ene%vdw_tot + &
                  pot_ene%elec_tot + &
                  pot_ene%gb + &
                  pot_ene%surf + &
                  pot_ene%bond + &
                  pot_ene%angle + &
                  pot_ene%dihedral + &
                  pot_ene%vdw_14 + &
                  pot_ene%elec_14 + &
                  pot_ene%restraint + &
                  pot_ene%angle_ub + &
                  pot_ene%imp + &
                  pot_ene%cmap


#ifdef MPI
! Distribute energies to all processes.

  call gb_distribute_enes(pot_ene)

#endif /* MPI */

 if(iamd.gt.1)then
 ! Calculate the boosting weight for dihedrals amd
   call calculate_amd_dih_weights(atm_cnt,pot_ene%dihedral,frc,crd)
 endif

#ifdef MPI
! Distribute forces to atom owners; this is a reduce sum operation:

  call gb_frcs_distrib(atm_cnt, frc, dbl_mpi_recv_buf)

#endif /* MPI */

  if (charmm_active .and. do_charmm_dump_gold == 1) then
#ifdef MPI
    call gb_mpi_gathervec(atm_cnt, frc)
    if (master) then
      call gb_charmm_dump_gold(atm_cnt,frc,pot_ene)
      write(mdout, '(a)') 'charmm_gold() completed. Exiting'
    endif
#else
    call gb_charmm_dump_gold(atm_cnt,frc,pot_ene)
    write(mdout, '(a)') 'charmm_gold() completed. Exiting'
#endif
    call mexit(6, 0)
  endif

! Calculate the NMR restraint energy contributions, if requested.
 
  if (nmropt .ne. 0) then
    call nmr_calc(crd, frc, enmr, 6)
    pot_ene%restraint = pot_ene%restraint + enmr(1) + enmr(2) + enmr(3)
    pot_ene%total = pot_ene%total + enmr(1) + enmr(2) + enmr(3)
  end if

!AMD DUAL BOOST CALC START
   if(iamd.gt.0)then
!calculate totboost and apply weight to frc. frc=frc*fwgt
     call calculate_amd_total_weights_gb(atm_cnt,pot_ene%total,pot_ene%dihedral,frc,crd)
   end if

! If belly is on then set the belly atom forces to zero:

  if (ibelly .gt. 0) call bellyf(atm_cnt, atm_igroup, frc)
#endif /* CUDA */
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
  use angles_ub_mod
  use bonds_mod
  use constraints_mod
  use dihedrals_mod
  use dihedrals_imp_mod
  use cmap_mod
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
  double precision              :: ub_ene
  double precision              :: dihedral_imp_ene
  double precision              :: cmap_ene

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

! UB Angle energy contribution
  if (cit_nthet_ub .gt. 0) then
    call get_angle_ub_energy(cit_nthet_ub, cit_angle_ub, crd, frc, ub_ene)
    pot_ene%angle_ub = ub_ene
  endif

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

! Improper dihedral energy contributions:

  if (cit_nimphi .gt. 0) then
    call get_dihed_imp_energy(cit_nimphi, cit_dihed_imp, crd, frc, dihedral_imp_ene)
    pot_ene%imp = dihedral_imp_ene
  endif

!CMAP
  if (cit_cmap_term_count .gt. 0) then
    call get_cmap_energy(cit_cmap_term_count, cit_cmap, crd, frc, cmap_ene)
    pot_ene%cmap = cmap_ene
  endif

  

! Dihedral energy contribution:

 if(iamd.gt.1)then
  if (ntf .le. 5) then
    if (cit_nphih .gt. 0) then
      call get_dihed_energy_amd(cit_nphih, cit_dihed, crd, dihedral_ene)
      pot_ene%dihedral = dihedral_ene
    end if
  end if

  if (ntf .le. 6) then
    if (cit_nphia .gt. 0) then
      call get_dihed_energy_amd(cit_nphia, cit_dihed(cit_nphih + 1), crd, &
                            dihedral_ene)
      pot_ene%dihedral = pot_ene%dihedral + dihedral_ene
    end if
  end if
 else
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
 endif

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
                     mpi_sum, pmemd_comm, err_code_mpi)

  pot_ene = dat_out

  return

end subroutine gb_distribute_enes
#endif

end module gb_force_mod
