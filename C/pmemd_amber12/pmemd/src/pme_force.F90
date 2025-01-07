#include "copyright.i"

!*******************************************************************************
!
! Module: pme_force_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module pme_force_mod

  use gbl_datatypes_mod
  use prmtop_dat_mod, only : numextra
  use extra_pnts_nb14_mod
#ifdef TIME_TEST
  use timers_mod
#endif

  implicit none

  ! Potential energies, with breakdown, from pme.  This is intended to be the
  ! external interface to potential energy data produced by this module, in
  ! particular by subroutine pme_force().

  type pme_pot_ene_rec
    sequence
    double precision    :: total
    double precision    :: vdw_tot      ! total of dir, recip
    double precision    :: vdw_dir
    double precision    :: vdw_recip
    double precision    :: elec_tot     ! total of dir, recip, nb_adjust, self
    double precision    :: elec_dir
    double precision    :: elec_recip
    double precision    :: elec_nb_adjust
    double precision    :: elec_self
    double precision    :: hbond
    double precision    :: bond
    double precision    :: angle
    double precision    :: dihedral
    double precision    :: vdw_14
    double precision    :: elec_14
    double precision    :: restraint
    double precision    :: angle_ub !CHARMM Urey-Bradley Energy
    double precision    :: imp      !CHARMM Improper Energy
    double precision    :: cmap     !CHARMM CMAP
  end type pme_pot_ene_rec

  integer, parameter    :: pme_pot_ene_rec_size = 19

  type(pme_pot_ene_rec), parameter      :: null_pme_pot_ene_rec = &
    pme_pot_ene_rec(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                    0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                    0.d0,0.d0,0.d0)
                         

  ! Virials, with breakdown, from pme.  This is intended to be the
  ! structure for storing virials; it currently is not exported but could be.

  type pme_virial_rec
    sequence
    double precision    :: molecular(3, 3)
    double precision    :: atomic(3, 3)
    double precision    :: elec_direct(3, 3)
    double precision    :: elec_nb_adjust(3, 3)
    double precision    :: elec_recip(3, 3)
    double precision    :: elec_recip_vdw_corr(3, 3)
    double precision    :: elec_recip_self(3, 3)
    double precision    :: elec_14(3, 3)
    double precision    :: ep_frame(3, 3)       ! for extra points...
    double precision    :: eedvir               ! used in Darden's error est.
  end type pme_virial_rec

  integer, parameter    :: pme_virial_rec_size = 82

  type(pme_virial_rec), parameter      :: null_pme_virial_rec = &
    pme_virial_rec(9*0.d0, 9*0.d0, 9*0.d0, &
                   9*0.d0, 9*0.d0, 9*0.d0, &
                   9*0.d0, 9*0.d0, 9*0.d0, &
                   0.d0)

  double precision, allocatable :: img_frc(:,:)
  double precision, allocatable :: nb_frc(:,:)
  double precision, allocatable :: ips_frc(:,:)
  integer, allocatable, save            :: gbl_nvdwcls(:)

  ! The following variables don't need to be broadcast:

  integer, save         :: irespa = 0

  ! Molecular virial correction factor; used internally.

  double precision, save, private       :: molvir_netfrc_corr(3, 3)

contains

!*******************************************************************************
!
! Subroutine:  alloc_pme_force_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_pme_force_mem(ntypes, num_ints, num_reals)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer                       :: ntypes

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed

  allocate(gbl_nvdwcls(ntypes), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(gbl_nvdwcls)

  gbl_nvdwcls(:) = 0

  return

end subroutine alloc_pme_force_mem

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_pme_force_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_pme_force_dat(ntypes)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer       :: ntypes

! Local variables:

  integer       :: num_ints, num_reals  ! returned values discarded

  ! Nothing to broadcast.  We just allocate storage in the non-master nodes.

  if (.not. master) then
    num_ints = 0
    num_reals = 0
    call alloc_pme_force_mem(ntypes, num_ints, num_reals)
  end if

  ! The allocated data is not initialized from the master node.

  return

end subroutine bcast_pme_force_dat

!*******************************************************************************
!
! Subroutine:  alloc_force_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_force_mem(atm_cnt,num_reals,ips)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_reals

  integer                       :: ips

! Local variables:

  integer               :: alloc_failed

  if (ips .gt. 0)then
   allocate(img_frc(3, atm_cnt), &
            nb_frc(3, atm_cnt), &
            ips_frc(3, atm_cnt), &
            stat = alloc_failed)

   if (alloc_failed .ne. 0) call setup_alloc_error
 
   num_reals = num_reals + size(img_frc) + size(nb_frc) + size(ips_frc)
  else
   allocate(img_frc(3, atm_cnt), &
            nb_frc(3, atm_cnt), &
            stat = alloc_failed)
 
   if (alloc_failed .ne. 0) call setup_alloc_error
 
   num_reals = num_reals + size(img_frc) + size(nb_frc) 
  endif

  return

end subroutine alloc_force_mem

!*******************************************************************************
!
! Subroutine:  dealloc_force_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine dealloc_force_mem(ips)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer                       :: ips

  if (ips .gt. 0)then
    deallocate(img_frc, nb_frc, ips_frc)
  else
    deallocate(img_frc, nb_frc)
  endif

  return

end subroutine dealloc_force_mem
#endif

!*******************************************************************************
!
! Subroutine:  pme_force
!
! Description: <TBS>
!              
!*******************************************************************************

#ifdef MPI
subroutine pme_force(atm_cnt, crd, frc, img_atm_map, atm_img_map, &
                     my_atm_lst, new_list, need_pot_enes, need_virials, &
                     pot_ene, virial, ekcmt, pme_err_est)

  use angles_mod
  use angles_ub_mod
  use bonds_mod
  use constraints_mod
  use dihedrals_mod
  use dihedrals_imp_mod
  use cmap_mod
  use dynamics_dat_mod
  use dynamics_mod
#ifdef DIRFRC_EFS
  use ene_frc_splines_mod
#endif /* DIRFRC_EFS */
  use nb_exclusions_mod
  use pme_direct_mod
  use pme_recip_dat_mod
  use pme_blk_recip_mod
  use pme_slab_recip_mod
  use loadbal_mod
  use img_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use nb_pairlist_mod
  use nmr_calls_mod
  use parallel_dat_mod
  use parallel_mod
  use pbc_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use runfiles_mod
  use timers_mod
  use charmm_mod, only : charmm_active
  use mdin_debugf_dat_mod, only : do_charmm_dump_gold
  use nbips_mod
  use amd_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt
  double precision              :: crd(3, atm_cnt)
  double precision              :: frc(3, atm_cnt)
  integer                       :: img_atm_map(atm_cnt)
  integer                       :: atm_img_map(atm_cnt)
  integer                       :: my_atm_lst(*)
  logical                       :: new_list
  logical                       :: need_pot_enes
  logical                       :: need_virials
  type(pme_pot_ene_rec)         :: pot_ene
  double precision              :: virial(3)            ! Only used for MD
  double precision              :: ekcmt(3)             ! Only used for MD
  double precision              :: pme_err_est          ! Only used for MD

! Local variables:

  type(pme_virial_rec)          :: vir
  double precision              :: enmr(3)
  double precision              :: vir_vs_ene
  integer                       :: atm_lst_idx
  integer                       :: alloc_failed
  integer                       :: i, j
  logical                       :: params_may_change
  logical                       :: onstep
  double precision              :: net_frcs(3)

  double precision              :: evdwex, eelex
  double precision              :: temp_tot_dih
  integer                       :: buf_size

  call zero_time()
  call zero_pme_time()
#ifdef CUDA
  if (nmropt .ne. 0) then
    call nmr_weight(atm_cnt, crd, 6)
    pot_ene = null_pme_pot_ene_rec
    params_may_change = .true.
  else
    params_may_change = .false.
  end if

  if (use_pme /= 0) then
    call self(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
               params_may_change)
    if (vdwmeth .eq. 1) then
        call vdw_correction(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                            params_may_change)
    end if
    call update_time(pme_misc_timer)
    if (need_pot_enes) then ! DO NOT TOUCH THIS, SEE RUNMD.F90 TO SET need_pot_enes
      if((iamd .eq. 2).or.(iamd .eq. 3))then
        totdih = 0.0
        if(num_amd_lag .eq. 0)then
          temp_tot_dih = 0.0
          call gpu_calculate_amd_dihedral_energy(totdih)
          buf_size = 1
          call mpi_allreduce(totdih, temp_tot_dih, &
                         buf_size, mpi_double_precision, &
                         mpi_sum, pmemd_comm, err_code_mpi)
          totdih = temp_tot_dih
        end if
        call gpu_calculate_amd_dihedral_weight(totdih)
      end if
      call gpu_pme_ene(ew_coeff, uc_volume, pot_ene, enmr, virial, ekcmt)
      call update_time(nonbond_time)
      if (need_virials) then
        vir%molecular(1,1) = virial(1)
        vir%molecular(2,2) = virial(2)
        vir%molecular(3,3) = virial(3)
      end if
      call dist_enes_virs_netfrcs(pot_ene, vir, ekcmt, net_frcs, &
                                  need_pot_enes, need_virials)   
      if (need_virials) then                                
        virial(1) = vir%molecular(1,1)
        virial(2) = vir%molecular(2,2)
        virial(3) = vir%molecular(3,3)
      end if
    else
      call gpu_pme_force(ew_coeff, uc_volume, virial, ekcmt)
      call update_time(nonbond_time)
        
      if (need_virials) then 
        vir%molecular(1,1) = virial(1)
        vir%molecular(2,2) = virial(2)
        vir%molecular(3,3) = virial(3) 
        call dist_enes_virs_netfrcs(pot_ene, vir, ekcmt, net_frcs, &
                                  need_pot_enes, need_virials)                                  
        virial(1) = vir%molecular(1,1)
        virial(2) = vir%molecular(2,2)
        virial(3) = vir%molecular(3,3)
      end if    
    end if
    call gpu_gather_pme_forces()  
  else
    call update_time(pme_misc_timer)
    call ipsupdate(ntb)
    if (need_pot_enes) then ! DO NOT TOUCH THIS, SEE RUNMD.F90 TO SET need_pot_enes
      if((iamd .eq. 2).or.(iamd .eq. 3))then
        totdih = 0.0
        if(num_amd_lag .eq. 0)then
          temp_tot_dih = 0.0
          call gpu_calculate_amd_dihedral_energy(totdih)
          buf_size = 1
          call mpi_allreduce(totdih, temp_tot_dih, &
                         buf_size, mpi_double_precision, &
                         mpi_sum, pmemd_comm, err_code_mpi)
          totdih = temp_tot_dih
        end if
        call gpu_calculate_amd_dihedral_weight(totdih)
      end if
      call gpu_ips_ene(uc_volume, pot_ene, enmr, virial, ekcmt)
      call update_time(nonbond_time)
      if (need_virials) then
        vir%molecular(1,1) = virial(1)
        vir%molecular(2,2) = virial(2)
        vir%molecular(3,3) = virial(3)
      end if
      call dist_enes_virs_netfrcs(pot_ene, vir, ekcmt, net_frcs, &
                                  need_pot_enes, need_virials)   
      if (need_virials) then                                
        virial(1) = vir%molecular(1,1)
        virial(2) = vir%molecular(2,2)
        virial(3) = vir%molecular(3,3)
      end if
    else
      call gpu_ips_force(uc_volume, virial, ekcmt)
      call update_time(nonbond_time)
        
      if (need_virials) then 
        vir%molecular(1,1) = virial(1)
        vir%molecular(2,2) = virial(2)
        vir%molecular(3,3) = virial(3) 
        call dist_enes_virs_netfrcs(pot_ene, vir, ekcmt, net_frcs, &
                                  need_pot_enes, need_virials)                                  
        virial(1) = vir%molecular(1,1)
        virial(2) = vir%molecular(2,2)
        virial(3) = vir%molecular(3,3)
      end if    
    end if
    call gpu_gather_ips_forces()  
  end if
  if (nmropt .ne. 0) then
    call nmr_calc(crd, frc, enmr, 6)
  end if
!AMD calculate weight and scale forces
  if(iamd.gt.0)then
    call gpu_calculate_and_apply_amd_weights(pot_ene%total, pot_ene%dihedral, &
                                             num_amd_lag)
  endif

  call update_time(fcve_dist_time)

#else      

  call do_load_balancing(new_list, atm_cnt)

! Do weight changes, if requested.

  if (nmropt .ne. 0) call nmr_weight(atm_cnt, crd, 6)

! If no ene / force calcs are to be done, clear everything and bag out now.

  if (ntf .eq. 8) then
    pot_ene = null_pme_pot_ene_rec
    virial(:) = 0.d0
    ekcmt(:) = 0.d0
    pme_err_est = 0.d0
    frc(:,:) = 0.d0
    return
  end if

! Calculate the non-bonded contributions:

! Direct part of ewald plus vdw, hbond, pairlist setup and image claiming:

  if (ntp .gt. 0) call fill_tranvec(gbl_tranvec)

  ! The following encapsulation (save_imgcrds) seems to be necessary to
  ! prevent an optimization bug with the SGI f90 compiler.  Sigh...

  if (new_list) then
    call pme_list(atm_cnt, crd, atm_nb_maskdata, atm_nb_mask)
    call start_loadbal_timer
    call save_imgcrds(gbl_img_crd, gbl_saved_imgcrd, &
                      gbl_used_img_cnt, gbl_used_img_lst)
  else
    call start_loadbal_timer
    call adjust_imgcrds(atm_cnt, gbl_img_crd, img_atm_map, &
                        gbl_used_img_cnt, gbl_used_img_lst, &
                        gbl_saved_imgcrd, crd, gbl_atm_saved_crd, ntp)
  end if

! Zero energies that are stack or call parameters:

  pot_ene = null_pme_pot_ene_rec
  vir = null_pme_virial_rec

  virial(:) = 0.d0
  ekcmt(:) = 0.d0
  pme_err_est = 0.d0

! Zero internal energies, virials, etc.

  enmr(:) = 0.d0
  vir_vs_ene = 0.d0

  net_frcs(:) = 0.d0
  molvir_netfrc_corr(:,:) = 0.d0

  do j = 1, gbl_used_img_cnt
    i = gbl_used_img_lst(j)
    img_frc(:, i) = 0.d0
  end do

  ! We also need to zero anything on the extra used atom list since nonbonded
  ! forces for it will get sent to the atom owner.  We don't calc any such
  ! forces for these atoms, but they are on the send atom list in order to
  ! get their coordinates updated.

  call zero_extra_used_atm_img_frcs(img_frc)

! Don't do recip if PME is not invoked. Don't do it this step unless
! mod(irespa,nrepsa) = 0

  onstep = mod(irespa, nrespa) .eq. 0

  params_may_change = (nmropt .ne. 0)

! Calculate the exclusions correction for ips if using it 
! Clear the force array and the atom-ordered nonbonded force array.
! We delay copying the nonbonded forces into the force
! array in order to be able to schedule i/o later, and batch it up.

  evdwex = 0.0
  eelex = 0.0

  if( ips > 0 ) then

    do atm_lst_idx = 1, my_atm_cnt
      i = my_atm_lst(atm_lst_idx)
      frc(:,i) = 0.d0
      nb_frc(:,i) = 0.d0
      ips_frc(:,i) = 0.d0
    end do
    call eexips(evdwex,eelex,ips_frc,crd, &
      gbl_img_qterm,gbl_ipairs,atm_nb_maskdata, &
      atm_nb_mask,img_atm_map,my_atm_cnt,gbl_tranvec,my_atm_lst)

    if (new_list) call get_img_frc_ips_distribution(atm_cnt)
    call distribute_img_ips_frcs(atm_cnt, ips_frc)
  else
    do atm_lst_idx = 1, my_atm_cnt
      i = my_atm_lst(atm_lst_idx)
      frc(:,i) = 0.d0
      nb_frc(:,i) = 0.d0
    end do
  endif

  if (use_pme /= 0 .and. onstep) then

    ! Self energy:

    ! The small amount of time used here gets lumped with the recip stuff...

    if (gbl_frc_ene_task_lst(1) .eq. mytaskid) then
      call self(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
                params_may_change)
    end if

    ! Reciprocal energy:

    ! We intentionally keep the load balance counter running through the
    ! fft dist transposes in the recip code; synchronization will mess up the
    ! times a bit, but when not all nodes are doing recip calcs, it will
    ! really help with load balancing.

    if (i_do_recip) then
      call update_pme_time(pme_misc_timer)
      call update_loadbal_timer(elapsed_100usec_atmuser)
      if (block_fft .eq. 0) then
        call do_slab_pmesh_kspace(img_frc, pot_ene%elec_recip, vir%elec_recip)
      else
        call do_blk_pmesh_kspace(img_frc, pot_ene%elec_recip, vir%elec_recip)
      end if
      call update_loadbal_timer(elapsed_100usec_recipfrc)
      if (nrespa .gt. 1) call respa_scale(atm_cnt, img_frc, nrespa)
    end if

! Long range dispersion contributions:

! Continuum method:

    if (vdwmeth .eq. 1 .and. gbl_frc_ene_task_lst(1) .eq. mytaskid) then
      call vdw_correction(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                          params_may_change)
    end if

  end if      ! respa .and. use_pme

  call update_loadbal_timer(elapsed_100usec_atmuser)

! Direct part of ewald plus vdw, hbond, force and energy calculations:

  call update_pme_time(pme_misc_timer)
  if(ips > 0)then
    call get_nb_ips_energy(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                     gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                     pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                     vir%eedvir, vir%elec_direct)
  else
    call get_nb_energy(img_frc, gbl_img_crd, gbl_img_qterm, gbl_eed_cub, &
                     gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                     pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                     vir%eedvir, vir%elec_direct)
  end if

  call update_pme_time(dir_frc_sum_timer)
  call update_loadbal_timer(elapsed_100usec_dirfrc)

  call update_pme_time(pme_misc_timer)

! Calculate 1-4 electrostatic energies, forces:

  if (charmm_active) then
    if (need_virials) then
      call get_nb14_energy(atm_qterm, crd, nb_frc, atm_iac, typ_ico, &
                           gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                           pot_ene%elec_14, pot_ene%vdw_14, &
                           vir%elec_14)
    else
      call get_nb14_energy(atm_qterm, crd, nb_frc, atm_iac, typ_ico, &
                         gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                         pot_ene%elec_14, pot_ene%vdw_14)
    end if
  else 
    if (need_virials) then
      call get_nb14_energy(atm_qterm, crd, nb_frc, atm_iac, typ_ico, &
                           gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                           pot_ene%elec_14, pot_ene%vdw_14, &
                           vir%elec_14)
    else
      call get_nb14_energy(atm_qterm, crd, nb_frc, atm_iac, typ_ico, &
                         gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                         pot_ene%elec_14, pot_ene%vdw_14)
    end if
  endif

  call update_pme_time(dir_frc_sum_timer)

! Adjust energies, forces for masked out pairs:

  if (use_pme /= 0) then
  call nb_adjust(my_atm_cnt, my_atm_lst, &
                 atm_qterm, crd, gbl_nb_adjust_pairlst, gbl_eed_cub, &
                 nb_frc, pot_ene%elec_nb_adjust, vir%elec_nb_adjust)
  endif

  call update_pme_time(adjust_masked_timer)

! Get COM-relative coordinates.
! Calculate total nonbonded energy components.

  pot_ene%vdw_tot = pot_ene%vdw_dir + &
                    pot_ene%vdw_recip

  pot_ene%elec_tot = pot_ene%elec_dir + &
                     pot_ene%elec_recip + &
                     pot_ene%elec_nb_adjust + &
                     pot_ene%elec_self

  if( ips > 0 ) then
   !add the IPS contribution
   pot_ene%vdw_tot = pot_ene%vdw_tot + evdwex
   pot_ene%elec_tot = pot_ene%elec_tot + eelex
  end if

  if (need_virials) then
    call get_atm_rel_crd(my_mol_cnt, gbl_mol_com, crd, atm_rel_crd, &
                         gbl_my_mol_lst)
  end if

  call update_pme_time(pme_misc_timer)
  call update_time(nonbond_time)

! Calculate the other contributions:

  call pme_bonded_force(crd, frc, pot_ene)

! The above stuff gets lumped as "other time"...

  ! Sum up total potential energy for this task:

  pot_ene%total = pot_ene%vdw_tot + &
                  pot_ene%elec_tot + &
                  pot_ene%hbond + &
                  pot_ene%bond + &
                  pot_ene%angle + &
                  pot_ene%dihedral + &
                  pot_ene%vdw_14 + &
                  pot_ene%elec_14 + &
                  pot_ene%restraint + &
                  pot_ene%imp + &
                  pot_ene%angle_ub + &
                  pot_ene%cmap
                               
  call update_loadbal_timer(elapsed_100usec_atmowner)

!add the ips non bonded contribution
  if( ips > 0 ) then
    do atm_lst_idx = 1, my_atm_cnt
      i = my_atm_lst(atm_lst_idx)
      nb_frc(:,i) = nb_frc(:,i) + ips_frc(:,i) 
    end do
  endif

  if (new_list) call get_img_frc_distribution(atm_cnt)

  call distribute_img_frcs(atm_cnt, img_frc, nb_frc)

  call update_time(fcve_dist_time)

  call zero_pme_time()

  ! If using extra points and a frame (checked internal to subroutine),
  ! transfer force and torque from the extra points to the parent atom:

  if (numextra .gt. 0 .and. frameon .ne. 0) &
    call orient_frc(crd, nb_frc, vir%ep_frame, ep_frames, gbl_frame_cnt)

  ! If the net force correction is in use, here we determine the net forces by
  ! looking at the sum of all nonbonded forces.  This should give the same
  ! result as just looking at the reciprocal forces, but it is more
  ! computationally convenient, especially for extra points, to do it this way.

  if (netfrc .gt. 0 .and. onstep) then
    do atm_lst_idx = 1, my_atm_cnt
      i = my_atm_lst(atm_lst_idx)
      net_frcs(:) = net_frcs(:) + nb_frc(:, i)
    end do
  end if

  ! First phase of virial work.  We need just the nonbonded forces at this
  ! stage.

  if (need_virials) then

    ! The relative coordinates for the extra points will mess up the molvir
    ! netfrc correction here unless we zero them out.  They are of no import
    ! in ekcom calc because extra points have a mass of 0.

    if (numextra .gt. 0 .and. frameon .ne. 0) &
      call zero_extra_pnts_vec(atm_rel_crd, ep_frames, gbl_frame_cnt)

    do atm_lst_idx = 1, my_atm_cnt
      i = my_atm_lst(atm_lst_idx)
      vir%molecular(:,1) = vir%molecular(:,1) + nb_frc(:,i) * atm_rel_crd(1,i)
      vir%molecular(:,2) = vir%molecular(:,2) + nb_frc(:,i) * atm_rel_crd(2,i)
      vir%molecular(:,3) = vir%molecular(:,3) + nb_frc(:,i) * atm_rel_crd(3,i)
      molvir_netfrc_corr(:, 1) = molvir_netfrc_corr(:, 1) + atm_rel_crd(1, i)
      molvir_netfrc_corr(:, 2) = molvir_netfrc_corr(:, 2) + atm_rel_crd(2, i)
      molvir_netfrc_corr(:, 3) = molvir_netfrc_corr(:, 3) + atm_rel_crd(3, i)
    end do

  ! Finish up virial work; Timing is inconsequential...

    vir%atomic(:,:) = vir%elec_recip(:,:) + &
                      vir%elec_direct(:,:) + &
                      vir%elec_nb_adjust(:,:) + &
                      vir%elec_recip_vdw_corr(:,:) + &
                      vir%elec_recip_self(:,:) + &
                      vir%elec_14(:,:) + &
                      vir%ep_frame(:,:)

!add the ips contribution to virials
    if( ips > 0 ) then
     !add the IPS contribution to virials
     vir%atomic(:,:) = vir%atomic(:,:) + VIREXIPS(:,:)
    endif

    vir%molecular(:,:) = vir%molecular(:,:) + vir%atomic(:,:)

    call get_ekcom(my_mol_cnt, gbl_mol_mass_inv, ekcmt, atm_vel, &
                   atm_mass, gbl_my_mol_lst)
  end if

  call update_time(nonbond_time)
  call update_pme_time(pme_misc_timer)

  ! Add a bunch of values (energies, virials, etc.) together as needed.

  call dist_enes_virs_netfrcs(pot_ene, vir, ekcmt, net_frcs, &
                              need_pot_enes, need_virials)
  
  call update_time(fcve_dist_time)
  call zero_pme_time()

  if(iamd.gt.1)then
    call calculate_amd_dih_weights(atm_cnt,pot_ene%dihedral,frc,crd)
  endif

  ! Copy image forces to atom forces, correcting for netfrc as you go, if
  ! appropriate.  Do not remove net force if netfrc = 0; e.g. in minimization.

  if (netfrc .gt. 0 .and. onstep) then

    net_frcs(:) = net_frcs(:) / dble(atm_cnt - numextra)

    do atm_lst_idx = 1, my_atm_cnt
      i = my_atm_lst(atm_lst_idx)
      frc(:, i) = frc(:, i) + nb_frc(:, i) - net_frcs(:)
    end do

    if (numextra .gt. 0 .and. frameon .ne. 0) &
      call zero_extra_pnts_vec(frc, ep_frames, gbl_frame_cnt)

    ! Correct the molecular virial for netfrc:

    if (need_virials) then
      vir%molecular(1,:) = vir%molecular(1,:) - &
                           molvir_netfrc_corr(1,:) * net_frcs(1)
      vir%molecular(2,:) = vir%molecular(2,:) - &
                           molvir_netfrc_corr(2,:) * net_frcs(2)
      vir%molecular(3,:) = vir%molecular(3,:) - &
                           molvir_netfrc_corr(3,:) * net_frcs(3)
    end if
  
  else

    do atm_lst_idx = 1, my_atm_cnt
      i = my_atm_lst(atm_lst_idx)
      frc(1, i) = frc(1, i) + nb_frc(1, i)
      frc(2, i) = frc(2, i) + nb_frc(2, i)
      frc(3, i) = frc(3, i) + nb_frc(3, i)
    end do
    
  end if

! Calculate vir_vs_ene in master.

  if (master) then

      vir_vs_ene = vir%elec_recip(1, 1) + &
                   vir%elec_recip(2, 2) + &
                   vir%elec_recip(3, 3) + &
                   vir%eedvir + &
                   vir%elec_nb_adjust(1, 1) + &
                   vir%elec_nb_adjust(2, 2) + &
                   vir%elec_nb_adjust(3, 3)

      ! Avoid divide-by-zero for pure neutral systems (l-j spheres), or if
      ! energies are not being calculated this step...

      if (pot_ene%elec_tot .ne. 0.0d0) then
        vir_vs_ene = abs(vir_vs_ene + pot_ene%elec_tot)/abs(pot_ene%elec_tot)
      else
        vir_vs_ene = 0.0d0
      end if
      pme_err_est = vir_vs_ene

  end if

  ! Save virials in form used in runmd:

  if (need_virials) then
    virial(1) = 0.5d0 * vir%molecular(1, 1)
    virial(2) = 0.5d0 * vir%molecular(2, 2)
    virial(3) = 0.5d0 * vir%molecular(3, 3)
  end if

! Calculate the NMR restraint energy contributions, if requested.
 
  if (nmropt .ne. 0) then
    call nmr_calc(crd, frc, enmr, 6)
    pot_ene%restraint = pot_ene%restraint + enmr(1) + enmr(2) + enmr(3)
    pot_ene%total = pot_ene%total + enmr(1) + enmr(2) + enmr(3)
  end if

  if (master .and. verbose .gt. 0) then
    call write_netfrc(net_frcs)
    call pme_verbose_print(pot_ene, vir, vir_vs_ene)
  end if

!AMD DUAL BOOST CALC START
   if(iamd.gt.0)then
!calculate totboost and apply weight to frc. frc=frc*fwgt
     call calculate_amd_total_weights(atm_cnt, pot_ene%total, &
                        pot_ene%dihedral, frc, crd, my_atm_lst)
   end if

! If belly is on then set the belly atom forces to zero:

  if (ibelly .gt. 0) call bellyf(atm_cnt, atm_igroup, frc)

  call update_time(nonbond_time)
  call update_pme_time(pme_misc_timer)

  if (charmm_active .and. do_charmm_dump_gold == 1) then
    call mpi_gathervec(atm_cnt, frc)
    if (master) then
      call pme_charmm_dump_gold(atm_cnt,frc,pot_ene)
      write(mdout, '(a)') 'charmm_gold() completed. Exiting'
    endif
    call mexit(6, 0)
  end if

#endif /* CUDA */

  return

end subroutine pme_force

#else /* uniprocessor version below: */

subroutine pme_force(atm_cnt, crd, frc, img_atm_map, atm_img_map, &
                     new_list, need_pot_enes, need_virials, &
                     pot_ene, virial, ekcmt, pme_err_est)

  use angles_mod
  use bonds_mod
  use constraints_mod
  use dihedrals_mod
  use dynamics_dat_mod
  use dynamics_mod
#ifdef DIRFRC_EFS
  use ene_frc_splines_mod
#endif /* DIRFRC_EFS */
  use nb_exclusions_mod
  use pme_direct_mod
  use pme_slab_recip_mod
  use loadbal_mod
  use img_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use nb_pairlist_mod
  use nmr_calls_mod
  use parallel_dat_mod
  use parallel_mod
  use pbc_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use runfiles_mod
  use timers_mod
  use charmm_mod, only : charmm_active
  use mdin_debugf_dat_mod, only : do_charmm_dump_gold
  use nbips_mod
  use amd_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt
  double precision              :: crd(3, atm_cnt)
  double precision              :: frc(3, atm_cnt)
  integer                       :: img_atm_map(atm_cnt)
  integer                       :: atm_img_map(atm_cnt)
  logical                       :: new_list
  logical                       :: need_pot_enes
  logical                       :: need_virials
  type(pme_pot_ene_rec)         :: pot_ene
  double precision, optional    :: virial(3)            ! Only used for MD
  double precision, optional    :: ekcmt(3)             ! Only used for MD
  double precision, optional    :: pme_err_est          ! Only used for MD

! Local variables:

  type(pme_virial_rec)          :: vir
  double precision              :: enmr(3)
  double precision              :: vir_vs_ene
  integer                       :: alloc_failed
  integer                       :: i, j
  logical                       :: params_may_change
  logical                       :: onstep
  double precision              :: net_frcs(3)
  double precision, allocatable :: img_frc(:,:)
  double precision              :: evdwex, eelex

  call zero_time()
  call zero_pme_time()
  onstep = mod(irespa, nrespa) .eq. 0
#ifdef CUDA
  if (nmropt .ne. 0) then
    call nmr_weight(atm_cnt, crd, 6)
    pot_ene = null_pme_pot_ene_rec
    params_may_change = .true.
    enmr(:) = 0.d0
  else
    params_may_change = .false.
  end if
  call update_time(pme_misc_timer)
  if (use_pme /= 0) then
    call self(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
               params_may_change)
    if (vdwmeth .eq. 1) then
      call vdw_correction(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                          params_may_change)
    end if
    if (need_pot_enes) then ! DO NOT TOUCH THIS, SEE RUNMD.F90 TO SET need_pot_enes
      if((iamd.eq.2).or.(iamd.eq.3))then
        call gpu_calculate_amd_dihedral_energy_weight()
      endif
      call gpu_pme_ene(ew_coeff, uc_volume, pot_ene, enmr, virial, ekcmt)
    else
      call gpu_pme_force(ew_coeff, uc_volume, virial, ekcmt)
    end if
  else
    call ipsupdate(ntb)
    if (need_pot_enes) then ! DO NOT TOUCH THIS, SEE RUNMD.F90 TO SET need_pot_enes
      if((iamd.eq.2).or.(iamd.eq.3))then
        call gpu_calculate_amd_dihedral_energy_weight()
      endif
      call gpu_ips_ene(uc_volume, pot_ene, enmr, virial, ekcmt)
    else
      call gpu_ips_force(uc_volume, virial, ekcmt) 
    end if
  end if
  if (nmropt .ne. 0) then
    call nmr_calc(crd, frc, enmr, 6)
  end if
  if(iamd.gt.0)then
!check if the energy here is updated
!AMD calculate weight and scale forces
    call gpu_calculate_and_apply_amd_weights(pot_ene%total, pot_ene%dihedral,&
                                             num_amd_lag)
  endif
  call update_time(nonbond_time)
#else    

! Zero energies that are stack or call parameters:

  pot_ene = null_pme_pot_ene_rec
  vir = null_pme_virial_rec

  virial(:) = 0.d0
  ekcmt(:) = 0.d0
  pme_err_est = 0.d0

! Zero internal energies, virials, etc.

  enmr(:) = 0.d0
  vir_vs_ene = 0.d0

  net_frcs(:) = 0.d0
  molvir_netfrc_corr(:,:) = 0.d0

! Do weight changes, if requested.

  if (nmropt .ne. 0) call nmr_weight(atm_cnt, crd, 6)

! If no force calcs are to be done, clear the frc array and bag out now.

  if (ntf .eq. 8) then
    frc(:,:) = 0.d0
    return
  end if

  allocate(img_frc(3, atm_cnt), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

! Calculate the non-bonded contributions:

! Direct part of ewald plus vdw, hbond, pairlist setup and image claiming:

  params_may_change = (nmropt .ne. 0)

  if (ntp .gt. 0) call fill_tranvec(gbl_tranvec)

  ! The following encapsulation (save_imgcrds) seems to be necessary to
  ! prevent an optimization bug with the SGI f90 compiler.  Sigh...

  if (new_list) then
    call pme_list(atm_cnt, crd, atm_nb_maskdata, atm_nb_mask)
    call save_imgcrds(atm_cnt, gbl_img_crd, gbl_saved_imgcrd)
  else
    call adjust_imgcrds(atm_cnt, gbl_img_crd, img_atm_map, &
                        gbl_saved_imgcrd, crd, gbl_atm_saved_crd, ntp)
  end if

  do i = 1, atm_cnt
    img_frc(1, i) = 0.d0
    img_frc(2, i) = 0.d0
    img_frc(3, i) = 0.d0
  end do
  
  evdwex = 0.0
  eelex = 0.0

  if( ips > 0 ) then
    frc(:,:) = 0.d0
    call eexips(evdwex,eelex,frc,crd, &
                gbl_img_qterm,gbl_ipairs,atm_nb_maskdata, &
                atm_nb_mask,img_atm_map,atm_cnt,gbl_tranvec)
  endif

! Don't do recip if PME is not invoked. Don't do it this step unless
! mod(irespa,nrepsa) = 0

  if (use_pme /= 0 .and. onstep) then

    ! Self energy:

    ! The small amount of time used here gets lumped with the recip stuff...

    call self(pot_ene%elec_self, ew_coeff, uc_volume, vir%elec_recip_self, &
              params_may_change)

    ! Reciprocal energy:

    call update_pme_time(pme_misc_timer)
    ! Only the slab fft implementation is used in uniprocessor runs...
    call do_slab_pmesh_kspace(img_frc, pot_ene%elec_recip, vir%elec_recip)
    if (nrespa .gt. 1) call respa_scale(atm_cnt, img_frc, nrespa)

! Long range dispersion contributions:

! Continuum method:

    if (vdwmeth .eq. 1) then
      call vdw_correction(pot_ene%vdw_recip, vir%elec_recip_vdw_corr, &
                          params_may_change)
    end if

  end if      ! respa

! Direct part of ewald plus vdw, hbond, force and energy calculations:

  call update_pme_time(pme_misc_timer)

  if(ips /= 0)then
    call get_nb_ips_energy(img_frc, gbl_img_crd, gbl_img_qterm,  gbl_eed_cub, &
                     gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                     pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                     vir%eedvir, vir%elec_direct)

  else
    call get_nb_energy(img_frc, gbl_img_crd, gbl_img_qterm,  gbl_eed_cub, &
                     gbl_ipairs, gbl_tranvec, need_pot_enes, need_virials, &
                     pot_ene%elec_dir, pot_ene%vdw_dir, pot_ene%hbond, &
                     vir%eedvir, vir%elec_direct)
  end if

  call update_pme_time(dir_frc_sum_timer)

! Transfer image forces to the force array.


  if( ips > 0 ) then
  do i = 1, atm_cnt
    j = atm_img_map(i)
    frc(1, i) = frc(1, i) + img_frc(1, j)
    frc(2, i) = frc(2, i) + img_frc(2, j)
    frc(3, i) = frc(3, i) + img_frc(3, j)
  end do
  else
  do i = 1, atm_cnt
    j = atm_img_map(i)
    frc(1, i) = img_frc(1, j)
    frc(2, i) = img_frc(2, j)
    frc(3, i) = img_frc(3, j)
  end do
  endif

  call update_pme_time(pme_misc_timer)

! Calculate 1-4 electrostatic energies, forces:

  if (charmm_active) then
    if (need_virials) then
      call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
                           gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                           pot_ene%elec_14, pot_ene%vdw_14, &
                           vir%elec_14)
    else
      call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
                           gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
                           pot_ene%elec_14, pot_ene%vdw_14)
    end if
  else
    if (need_virials) then
      call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
                           gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                           pot_ene%elec_14, pot_ene%vdw_14, &
                           vir%elec_14)
    else
      call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
                           gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
                           pot_ene%elec_14, pot_ene%vdw_14)
    end if
  end if

  call update_pme_time(dir_frc_sum_timer)

! Adjust energies, forces for masked out pairs:

  if (use_pme /= 0) then
  call nb_adjust(atm_cnt, atm_qterm, crd, gbl_nb_adjust_pairlst, &
                 gbl_eed_cub, frc, pot_ene%elec_nb_adjust, vir%elec_nb_adjust)
  endif

  call update_pme_time(adjust_masked_timer)

  ! If using extra points and a frame (checked internal to subroutine),
  ! transfer force and torque from the extra points to the parent atom:

  if (numextra .gt. 0 .and. frameon .ne. 0) &
    call orient_frc(crd, frc, vir%ep_frame, ep_frames, gbl_frame_cnt)

  call update_pme_time(dir_frc_sum_timer)

  ! Calculate total nonbonded energy components.

  pot_ene%vdw_tot = pot_ene%vdw_dir + &
                    pot_ene%vdw_recip

  pot_ene%elec_tot = pot_ene%elec_dir + &
                     pot_ene%elec_recip + &
                     pot_ene%elec_nb_adjust + &
                     pot_ene%elec_self

  if( ips > 0 ) then
   !add the IPS contribution
   pot_ene%vdw_tot = pot_ene%vdw_tot + evdwex
   pot_ene%elec_tot = pot_ene%elec_tot + eelex
  end if

  ! If the net force correction is in use, here we determine the net forces by
  ! looking at the sum of all nonbonded forces.  This should give the same
  ! result as just looking at the reciprocal forces, but it is more
  ! computationally convenient, especially for extra points, to do it this way.

  ! NOTE: this is not currently done for the GPU

  if (netfrc .gt. 0 .and. onstep) then

    do i = 1, atm_cnt
      net_frcs(:) = net_frcs(:) + frc(:, i)
    end do

    ! Now do the correction:

    net_frcs(:) = net_frcs(:) / dble(atm_cnt - numextra)

    do i = 1, atm_cnt
      frc(:, i) = frc(:, i) - net_frcs(:)
    end do

    ! Any extra points must have their 0.d0 forces reset...

    if (numextra .gt. 0 .and. frameon .ne. 0) &
      call zero_extra_pnts_vec(frc, ep_frames, gbl_frame_cnt)

  end if

  call update_pme_time(pme_misc_timer)

  ! First phase of virial work.  We need just the nonbonded forces at this
  ! stage, though there will be a further correction in get_dihed.

  if (need_virials) then

    call get_atm_rel_crd(my_mol_cnt, gbl_mol_com, crd, atm_rel_crd)
    do i = 1, atm_cnt
      vir%molecular(:,1) = vir%molecular(:,1) + frc(:,i) * atm_rel_crd(1,i)
      vir%molecular(:,2) = vir%molecular(:,2) + frc(:,i) * atm_rel_crd(2,i)
      vir%molecular(:,3) = vir%molecular(:,3) + frc(:,i) * atm_rel_crd(3,i)
    end do

    call get_ekcom(my_mol_cnt, gbl_mol_mass_inv, ekcmt, atm_vel, atm_mass)

  end if

  call update_time(nonbond_time)
  call update_pme_time(pme_misc_timer)

  ! Calculate the other contributions:

  call pme_bonded_force(crd, frc, pot_ene)

  ! Sum up total potential energy for this task:

  pot_ene%total = pot_ene%vdw_tot + &
                  pot_ene%elec_tot + &
                  pot_ene%hbond + &
                  pot_ene%bond + &
                  pot_ene%angle + &
                  pot_ene%dihedral + &
                  pot_ene%vdw_14 + &
                  pot_ene%elec_14 + &
                  pot_ene%restraint + &
                  pot_ene%imp + &
                  pot_ene%angle_ub + &
                  pot_ene%cmap

 if(iamd.gt.1)then

! Calculate the boosting weight for amd
   call calculate_amd_dih_weights(atm_cnt,pot_ene%dihedral,frc,crd)
 endif

  ! Adjustment of total energy for constraint energies does not seem
  ! consistent, but it matches sander...

  call zero_pme_time()
  
  ! Finish up virial work; Timing is inconsequential...

  if (need_virials) then

    vir%atomic(:,:) = vir%elec_recip(:,:) + &
                      vir%elec_direct(:,:) + &
                      vir%elec_nb_adjust(:,:) + &
                      vir%elec_recip_vdw_corr(:,:) + &
                      vir%elec_recip_self(:,:) + &
                      vir%elec_14(:,:) + &
                      vir%ep_frame(:,:)

    if( ips > 0 ) then
     !add the IPS contribution to virials
      vir%atomic(:,:) = vir%atomic(:,:) + VIREXIPS(:,:)
    endif

    vir%molecular(:,:) = vir%molecular(:,:) + vir%atomic(:,:)

    ! Save virials in form used in runmd:

    virial(1) = 0.5d0 * vir%molecular(1, 1)
    virial(2) = 0.5d0 * vir%molecular(2, 2)
    virial(3) = 0.5d0 * vir%molecular(3, 3)

  end if

  vir_vs_ene = vir%elec_recip(1, 1) + &
               vir%elec_recip(2, 2) + &
               vir%elec_recip(3, 3) + &
               vir%eedvir + &
               vir%elec_nb_adjust(1, 1) + &
               vir%elec_nb_adjust(2, 2) + &
               vir%elec_nb_adjust(3, 3)

  ! Avoid divide-by-zero for pure neutral systems (l-j spheres), or if
  ! energies were not calculated this step...

  if (pot_ene%elec_tot .ne. 0.0d0) then
    vir_vs_ene = abs(vir_vs_ene + pot_ene%elec_tot)/abs(pot_ene%elec_tot)
  else
    vir_vs_ene = 0.0d0
  end if

  pme_err_est = vir_vs_ene

! Calculate the NMR restraint energy contributions, if requested.
 
  if (nmropt .ne. 0) then
    call nmr_calc(crd, frc, enmr, 6)
    pot_ene%restraint = pot_ene%restraint + enmr(1) + enmr(2) + enmr(3)
    pot_ene%total = pot_ene%total + enmr(1) + enmr(2) + enmr(3)
  end if

  if (verbose .gt. 0) then
    call write_netfrc(net_frcs)
    call pme_verbose_print(pot_ene, vir, vir_vs_ene)
  end if

!AMD DUAL BOOST CALC START
   if(iamd.gt.0)then
     call calculate_amd_total_weights(atm_cnt,pot_ene%total,pot_ene%dihedral,frc,crd)
   end if

! If belly is on then set the belly atom forces to zero:

  if (ibelly .gt. 0) call bellyf(atm_cnt, atm_igroup, frc)

  call update_time(nonbond_time)
  call update_pme_time(pme_misc_timer)
  if (charmm_active .and. do_charmm_dump_gold == 1) then
    call pme_charmm_dump_gold(atm_cnt,frc,pot_ene)
    write(mdout, '(a)') 'charmm_gold() completed. Exiting'
    call mexit(6, 0)
  end if

  deallocate(img_frc)
#endif /* CUDA */
  return

end subroutine pme_force

#endif /* MPI */

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  dist_enes_virs_netfrcs
!
! Description: We reduce the appropriate subset of values in the ene array,
!              the ekcmt array, and the pme_ene_vir common block.
!*******************************************************************************

subroutine dist_enes_virs_netfrcs(pot_ene, vir, ekcmt, net_frcs, &
                                  need_pot_enes, need_virials)

  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  type(pme_pot_ene_rec) :: pot_ene
  type(pme_virial_rec)  :: vir
  double precision      :: ekcmt(3)
  double precision      :: net_frcs(3)
  logical, intent(in)   :: need_pot_enes
  logical, intent(in)   :: need_virials

! Local variables:

  type pme_dat
    sequence
    type(pme_pot_ene_rec)       :: pot_ene
    type(pme_virial_rec)        :: vir
    double precision            :: ekcmt(3)
    double precision            :: molvir_netfrc_corr(3,3)
    double precision            :: net_frcs(3)
  end type pme_dat

  integer, parameter            :: pme_dat_size = &
                                   pme_pot_ene_rec_size + &
                                   pme_virial_rec_size + 3 + 9 + 3

  integer, parameter            :: pme_dat_no_pot_ene_size = &
                                   pme_virial_rec_size + 3 + 9 + 3
 
  type(pme_dat), save           :: dat_in               ! used by all tasks
  type(pme_dat), save           :: dat_out              ! used by master

  integer                       :: buf_size

#ifdef COMM_TIME_TEST
  call start_test_timer(6, 'dist_enes_virs_netfrcs', 0)
#endif

  if (need_pot_enes .or. nmropt .ne. 0 .or. verbose .ne. 0) then

    ! When you need potential energies, you need parts of the virial info in
    ! the master and adding in net frc info is then a minor impact, so send
    ! it all.

    dat_in%pot_ene                 = pot_ene
    dat_in%vir                     = vir
    dat_in%ekcmt(:)                = ekcmt(:)
    dat_in%molvir_netfrc_corr(:,:) = molvir_netfrc_corr(:,:)
    dat_in%net_frcs(:)             = net_frcs(:)

    buf_size = pme_dat_size

    call mpi_allreduce(dat_in%pot_ene%total, dat_out%pot_ene%total, &
                       buf_size, mpi_double_precision, &
                       mpi_sum, pmemd_comm, err_code_mpi)

    pot_ene                           = dat_out%pot_ene
    vir                               = dat_out%vir
    ekcmt(:)                          = dat_out%ekcmt(:)
    molvir_netfrc_corr(:,:)           = dat_out%molvir_netfrc_corr(:,:)
    net_frcs(:)                       = dat_out%net_frcs(:)

  else if (need_virials) then

    ! The virials are distributed, with net frc info added as it has
    ! a minor impact.

    dat_in%vir                     = vir
    dat_in%ekcmt(:)                = ekcmt(:)
    dat_in%molvir_netfrc_corr(:,:) = molvir_netfrc_corr(:,:)
    dat_in%net_frcs(:)             = net_frcs(:)

    buf_size = pme_dat_no_pot_ene_size

    call mpi_allreduce(dat_in%vir%molecular(1, 1), &
                       dat_out%vir%molecular(1, 1), &
                       buf_size, mpi_double_precision, &
                       mpi_sum, pmemd_comm, err_code_mpi)

    vir                               = dat_out%vir
    ekcmt(:)                          = dat_out%ekcmt(:)
    molvir_netfrc_corr(:,:)           = dat_out%molvir_netfrc_corr(:,:)
    net_frcs(:)                       = dat_out%net_frcs(:)

  else if (netfrc .ne. 0) then

    ! We just need net frc info!

    dat_in%net_frcs(:) = net_frcs(:)
    buf_size = 3

    call mpi_allreduce(dat_in%net_frcs(1), &
                       dat_out%net_frcs(1), &
                       buf_size, mpi_double_precision, &
                       mpi_sum, pmemd_comm, err_code_mpi)

    net_frcs(:)                       = dat_out%net_frcs(:)

  end if

#ifdef COMM_TIME_TEST
  call stop_test_timer(6)
#endif

  return

end subroutine dist_enes_virs_netfrcs
#endif

!*******************************************************************************
!
! Subroutine:  self
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine self(ene, ewaldcof, vol, vir, params_may_change)

  use gbl_constants_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  double precision      :: ene, ewaldcof, vol, vir(3, 3)
  logical               :: params_may_change

! Local variables:

  integer                       :: i
  double precision              :: ee_plasma
  logical, save                 :: setup_not_done = .true.
  double precision, save        :: factor
  double precision, save        :: sqrt_pi
  double precision, save        :: sumq
  double precision, save        :: sumq2

! Only compute sumq and sumq2 at beginning. They don't change. This code is
! only executed by the master, so we precalc anything we can...

  if (setup_not_done .or. params_may_change) then

    factor = -0.5d0 * PI / (ewaldcof * ewaldcof)

    sqrt_pi = sqrt(PI)

    sumq = 0.d0
    sumq2 = 0.d0

    do i = 1, natom
      sumq = sumq + atm_qterm(i)
      sumq2 = sumq2 + atm_qterm(i) * atm_qterm(i)
    end do

    setup_not_done = .false.

  end if

  ee_plasma = factor * sumq * sumq / vol

  ene = - sumq2 * ewaldcof / sqrt_pi + ee_plasma

#ifdef CUDA
  call gpu_self(ee_plasma, ene)
#endif

  ! The off-diagonal elements are already zero.

  vir(1,1) = -ee_plasma
  vir(2,2) = -ee_plasma
  vir(3,3) = -ee_plasma

  return

end subroutine self

!*******************************************************************************
!
! Subroutine:  vdw_correction
!
! Description:  Get analytic estimate of energy and virial corrections due to
!               dispersion interactions beyond the cutoff.
!*******************************************************************************

subroutine vdw_correction(ene, virial, params_may_change)

  use mdin_ctrl_dat_mod
  use pbc_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  double precision      :: ene, virial(3, 3)
  logical               :: params_may_change

! Local variables:

  integer                       :: i, j, ic, iaci
  logical, save                 :: setup_not_done = .true.
  double precision, save        :: ene_factor            ! Result of precalc.
  double precision              :: prefac, term

! Only compute ene_factor at the beginning. It doesn't change. This code is
! only executed by the master, so we precalc anything we can...

  if (setup_not_done .or. params_may_change) then

    term = 0.d0

    ! Will later divide by volume, which is all that could change:

    prefac = 2.d0 * PI / (3.d0 * vdw_cutoff**3)

    do i = 1, ntypes
      iaci = ntypes * (i - 1)
      do j = 1, ntypes
        ic = typ_ico(iaci + j)
        if (ic .gt. 0) term = term + &
                              gbl_nvdwcls(i) * gbl_nvdwcls(j) * gbl_cn2(ic)
      end do
    end do

    ene_factor = -prefac * term

    setup_not_done = .false.

#ifdef CUDA
    call gpu_vdw_correction(ene_factor / uc_volume)
#endif

  end if

  ene = ene_factor / uc_volume

  ! The off-diagonal elements are already zero.

  virial(1, 1) = - 2.d0 * ene
  virial(2, 2) = - 2.d0 * ene
  virial(3, 3) = - 2.d0 * ene

  return

end subroutine vdw_correction

!*******************************************************************************
!
! Subroutine:  respa_scale
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine respa_scale(atm_cnt, img_frc, nrespa)

  use img_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: img_frc(3, atm_cnt)
  integer               :: nrespa

! Local variables:

  integer               :: img_idx

#ifdef MPI
  integer               :: lst_idx
#endif /* MPI */

#ifdef MPI
  do lst_idx = 1, gbl_used_img_cnt
    img_idx = gbl_used_img_lst(lst_idx)
#else
  do img_idx = 1, atm_cnt
#endif
    img_frc(1, img_idx) = nrespa * img_frc(1, img_idx)
    img_frc(2, img_idx) = nrespa * img_frc(2, img_idx)
    img_frc(3, img_idx) = nrespa * img_frc(3, img_idx)
  end do

  return

end subroutine respa_scale

!*******************************************************************************
!
! Subroutine:  pme_bonded_force
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine pme_bonded_force(crd, frc, pot_ene)

  use angles_mod
  use angles_ub_mod
  use bonds_mod
  use constraints_mod
  use dihedrals_mod
  use dihedrals_imp_mod
  use cmap_mod
  use dynamics_mod
  use dynamics_dat_mod
  use mdin_ctrl_dat_mod
  use nmr_calls_mod
  use parallel_dat_mod
  use timers_mod

  implicit none

! Formal arguments:

  double precision              :: crd(3, *)
  double precision              :: frc(3, *)
  type(pme_pot_ene_rec)         :: pot_ene

  ! These energy variables are temporaries, for summing. DON'T use otherwise!

  double precision              :: bond_ene
  double precision              :: ub_ene
  double precision              :: angle_ene
  double precision              :: dihedral_ene
  double precision              :: dihedral_imp_ene
  double precision              :: cmap_ene

#ifdef MPI
  if (my_atm_cnt .eq. 0) return
#endif /* MPI */

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

end subroutine pme_bonded_force


!*******************************************************************************
!
! Subroutine:  write_netfrc
!
! Description:  Get the netfrc's back into the external axis order and print
!               them out.  We do this all in a separate subroutine just to
!               keep from cluttering up critical code.
!              
!*******************************************************************************

subroutine write_netfrc(net_frcs)

  use axis_optimize_mod
  use file_io_dat_mod

  implicit none

! Formal arguments:

  double precision      :: net_frcs(3)

! Local variables:

  integer               :: ord1, ord2, ord3

  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

  write(mdout, 33) net_frcs(ord1), net_frcs(ord2), net_frcs(ord3)

  return

33     format(1x, 'NET FORCE PER ATOM: ', 3(1x, e12.4))

end subroutine write_netfrc

!*******************************************************************************
!
! Subroutine:  pme_verbose_print
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine pme_verbose_print(pot_ene, vir, vir_vs_ene)

  use axis_optimize_mod
  use file_io_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  type(pme_pot_ene_rec) :: pot_ene
  type(pme_virial_rec)  :: vir
  double precision      :: vir_vs_ene

! Local variables:

  integer       :: ord1, ord2, ord3

  if (.not. master) return

  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

  if (verbose .ge. 1) then
    write(mdout, '(4(/,5x,a,f22.12))') &
          'Evdw                   = ', pot_ene%vdw_tot, &
          'Ehbond                 = ', pot_ene%hbond, &
          'Ecoulomb               = ', pot_ene%elec_tot
    write(mdout, '(2(/,5x,a,f22.12))') &
          'Iso virial             = ',  &
          vir%molecular(1, 1) + vir%molecular(2, 2) + vir%molecular(3, 3), &
          'Eevir vs. Ecoulomb     = ', vir_vs_ene
  end if

  if (verbose .ge. 2) then
    write(mdout, '(4(/,5x,a,f22.12),/)') &
          'E electrostatic (self) = ', pot_ene%elec_self, &
          '                (rec)  = ', pot_ene%elec_recip, &
          '                (dir)  = ', pot_ene%elec_dir, &
          '                (adj)  = ', pot_ene%elec_nb_adjust
    write(mdout, 30) vir%molecular(ord1, ord1), &
                     vir%molecular(ord1, ord2), &
                     vir%molecular(ord1, ord3)
    write(mdout, 30) vir%molecular(ord2, ord1), &
                     vir%molecular(ord2, ord2), &
                     vir%molecular(ord2, ord3)
    write(mdout, 30) vir%molecular(ord3, ord1), &
                     vir%molecular(ord3, ord2), &
                     vir%molecular(ord3, ord3)
30     format(5x, 'MOLECULAR VIRIAL: ', 3(1x, e14.8))
  end if

  if (verbose .eq. 3) then
    write(mdout, *) '--------------------------------------------'
    write(mdout, 31) vir%elec_recip(ord1, ord1), &
                     vir%elec_recip(ord1, ord2), &
                     vir%elec_recip(ord1, ord3)
    write(mdout, 31) vir%elec_recip(ord2, ord1), &
                     vir%elec_recip(ord2, ord2), &
                     vir%elec_recip(ord2, ord3)
    write(mdout, 31) vir%elec_recip(ord3, ord1), &
                     vir%elec_recip(ord3, ord2), &
                     vir%elec_recip(ord3, ord3)
    write(mdout, *) '..................'
31     format(5x, 'Reciprocal VIRIAL: ', 3(1x, e14.8))
    write(mdout, 32) vir%elec_direct(ord1, ord1), &
                     vir%elec_direct(ord1, ord2), &
                     vir%elec_direct(ord1, ord3)
    write(mdout, 32) vir%elec_direct(ord2, ord1), &
                     vir%elec_direct(ord2, ord2), &
                     vir%elec_direct(ord2, ord3)
    write(mdout, 32) vir%elec_direct(ord3, ord1), &
                     vir%elec_direct(ord3, ord2), &
                     vir%elec_direct(ord3, ord3)
    write(mdout, *) '..................'
32     format(5x, 'Direct VIRIAL: ', 3(1x, e14.8))
    write(mdout, 38) vir%eedvir
    write(mdout, *) '..................'
38     format(5x, 'Dir Sum EE vir trace: ', e14.8)
    write(mdout, 33) vir%elec_nb_adjust(ord1, ord1), &
                     vir%elec_nb_adjust(ord1, ord2), &
                     vir%elec_nb_adjust(ord1, ord3)
    write(mdout, 33) vir%elec_nb_adjust(ord2, ord1), &
                     vir%elec_nb_adjust(ord2, ord2), &
                     vir%elec_nb_adjust(ord2, ord3)
    write(mdout, 33) vir%elec_nb_adjust(ord3, ord1), &
                     vir%elec_nb_adjust(ord3, ord2), &
                     vir%elec_nb_adjust(ord3, ord3)
    write(mdout, *) '..................'
33     format(5x, 'Adjust VIRIAL: ', 3(1x, e14.8))
    write(mdout, 34) vir%elec_recip_vdw_corr(ord1, ord1), &
                     vir%elec_recip_vdw_corr(ord1, ord2), &
                     vir%elec_recip_vdw_corr(ord1, ord3)
    write(mdout, 34) vir%elec_recip_vdw_corr(ord2, ord1), &
                     vir%elec_recip_vdw_corr(ord2, ord2), &
                     vir%elec_recip_vdw_corr(ord2, ord3)
    write(mdout, 34) vir%elec_recip_vdw_corr(ord3, ord1), &
                     vir%elec_recip_vdw_corr(ord3, ord2), &
                     vir%elec_recip_vdw_corr(ord3, ord3)
    write(mdout, *) '..................'
34     format(5x, 'Recip Disp. VIRIAL: ', 3(1x, e14.8))
    write(mdout, 35) vir%elec_recip_self(ord1, ord1), &
                     vir%elec_recip_self(ord1, ord2), &
                     vir%elec_recip_self(ord1, ord3)
    write(mdout, 35) vir%elec_recip_self(ord2, ord1), &
                     vir%elec_recip_self(ord2, ord2), &
                     vir%elec_recip_self(ord2, ord3)
    write(mdout, 35) vir%elec_recip_self(ord3, ord1), &
                     vir%elec_recip_self(ord3, ord2), &
                     vir%elec_recip_self(ord3, ord3)
    write(mdout, *) '..................'
35     format(5x, 'Self VIRIAL: ', 3(1x, e14.8))
    write(mdout, 36) vir%elec_14(ord1, ord1), &
                     vir%elec_14(ord1, ord2), &
                     vir%elec_14(ord1, ord3)
    write(mdout, 36) vir%elec_14(ord2, ord1), &
                     vir%elec_14(ord2, ord2), &
                     vir%elec_14(ord2, ord3)
    write(mdout, 36) vir%elec_14(ord3, ord1), &
                     vir%elec_14(ord3, ord2), &
                     vir%elec_14(ord3, ord3)
    write(mdout, *) '..................'
36     format(5x, 'E14 VIRIAL: ', 3(1x, e14.8))
    write(mdout, 37) vir%atomic(ord1, ord1), &
                     vir%atomic(ord1, ord2), &
                     vir%atomic(ord1, ord3)
    write(mdout, 37) vir%atomic(ord2, ord1), &
                     vir%atomic(ord2, ord2), &
                     vir%atomic(ord2, ord3)
    write(mdout, 37) vir%atomic(ord3, ord1), &
                     vir%atomic(ord3, ord2), &
                     vir%atomic(ord3, ord3)
37     format(5x, 'Atomic VIRIAL: ', 3(1x, e14.8))
    write(mdout, *)'--------------------------------------------'
  end if

  return

end subroutine pme_verbose_print

end module pme_force_mod
