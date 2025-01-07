#include "copyright.i"

!*******************************************************************************
!
! Module:  pme_direct_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module pme_direct_mod

  implicit none

  ! The following storage is per-process common; ie., it SHOULD be
  ! broadcast from the master to the other processes!

  integer, parameter    :: pme_direct_int_cnt = 1

  integer                  mxeedtab

  common / pme_direct_int / mxeedtab

  save  :: / pme_direct_int /

  double precision, allocatable, save   :: gbl_eed_cub(:)

contains

!*******************************************************************************
!
! Subroutine:  init_pme_direct_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine init_pme_direct_dat(num_ints, num_reals)

  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use prmtop_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  double precision      :: dxdr
  integer               :: i

! Setup dxdr map from r to x in table lookup of eed

  dxdr = ew_coeff
  
! For eed table, assume all nonbond distances are less than 1.5 * es_cutoff
! between nonbond updates; i.e. no excess motion; this is enforced by
! es_cutoff.

  mxeedtab = int(dxdr * eedtbdns * es_cutoff * 1.5d0)

  call alloc_pme_direct_mem(num_ints, num_reals)

  return

end subroutine init_pme_direct_dat

!*******************************************************************************
!
! Subroutine:  alloc_pme_direct_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_pme_direct_mem(num_ints, num_reals)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed

  allocate(gbl_eed_cub(4 * mxeedtab), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(gbl_eed_cub)

  gbl_eed_cub(:) = 0.d0

  return

end subroutine alloc_pme_direct_mem

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_pme_direct_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_pme_direct_dat

  use parallel_dat_mod

  implicit none

! Local variables:

  integer       :: num_ints, num_reals  ! returned values discarded

  call mpi_bcast(mxeedtab, pme_direct_int_cnt, mpi_integer, 0, &
                 mpi_comm_world, err_code_mpi)

  if (.not. master) then
    num_ints = 0
    num_reals = 0
    call alloc_pme_direct_mem(num_ints, num_reals)
  end if

  ! The allocated data is not initialized from the master node.

  return

end subroutine bcast_pme_direct_dat
#endif

!*******************************************************************************
!
! Subroutine:  pme_list
!
! Description:  Handles set-up and error checking for calling of get_nb_list
!               which creates the nonbond list.
!
!*******************************************************************************

subroutine pme_list(atm_cnt, crd, atm_maskdata, atm_mask)

  use cit_mod
  use constraints_mod
  use pmemd_lib_mod
  use pbc_mod
  use nb_pairlist_mod
  use img_mod
  use pme_blk_recip_mod
  use pme_slab_recip_mod
  use pme_recip_dat_mod
  use loadbal_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_dat_mod
  use parallel_mod
  use pbc_mod
  use prmtop_dat_mod
  use timers_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  type(listdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)

! Local variables:

  integer               :: ifail
  integer               :: alloc_failed
  logical               :: dont_skip_belly_pairs

  ! We 0-base the following array for efficiency, but don't use atm_lst(0) so
  ! we can use 0 as the distinguished non-entry (NULL).

  double precision      :: fraction(3, atm_cnt) ! in range 0.0 - +0.999...
  type(cit_tbl_rec)     :: crd_idx_tbl(0 : cit_tbl_x_dim - 1, &
                                       0 : cit_tbl_y_dim - 1, &
                                       0 : cit_tbl_z_dim - 1)

  ! Under amber 7, when both members of an atom pair are in the belly (ie.,
  ! fixed), we don't add them to the nonbonded pairlist.  This saves time but
  ! makes the energies look different (in a manner that actually does not
  ! affect the simulation).

  dont_skip_belly_pairs = ibelly .eq. 0
   
  gbl_saved_box(:) = pbc_box(:)     ! Needed when pressure scaling.

#ifndef MPI
  call save_all_atom_crds(atm_cnt, crd, gbl_atm_saved_crd)
#endif

  call get_fract_crds(atm_cnt, crd, fraction)

  call setup_cit(atm_cnt, fraction, crd_idx_tbl, &
                 gbl_atm_img_map, gbl_img_atm_map)
#ifdef MPI
  call init_used_img_map(gbl_used_img_map, gbl_used_img_cnt, gbl_used_img_lst)
#endif /* MPI */

  call update_pme_time(cit_setup_timer)

  ! NOTE - It is important to note that no image mapping should occur prior
  !        to running get_nb_list().  We now exit setup_cit() without the
  !        used image map set also.  This is all being done to facilitate
  !        speed in image mapping in get_nb_list(), which should always be
  !        called immediately after setup_cit().

#ifdef MPI
  call map_pairlist_imgs(atm_cnt, crd_idx_tbl, &
                         fraction, atm_qterm, atm_iac, gbl_img_crd, &
                         gbl_img_qterm, gbl_img_iac, gbl_mapped_img_cnt, &
                         gbl_mapped_img_lst, gbl_img_atm_map)
#else
  call map_pairlist_imgs(atm_cnt, fraction, atm_qterm, atm_iac, gbl_img_crd, &
                         gbl_img_qterm, gbl_img_iac, gbl_img_atm_map)
#endif

  do

    ifail = 0

    call get_nb_list(atm_cnt, crd_idx_tbl, gbl_img_crd, gbl_img_atm_map, &
#ifdef MPI
                     gbl_used_img_map, gbl_used_img_cnt, gbl_used_img_lst, &
#endif
                     typ_ico, fraction, gbl_tranvec, atm_maskdata, &
                     atm_mask, gbl_atm_img_map, gbl_excl_img_flags, &
                     gbl_img_iac, atm_igroup, ntypes, ibelly, &
                     es_cutoff, vdw_cutoff, skinnb, verbose, ifail)

    if (ifail .eq. 0) exit

    ! Deallocate the old array, grow it by 10%, reallocate,
    ! and go back up to the top to try again.

    deallocate(gbl_ipairs)

    ipairs_maxsize = 11 * ipairs_maxsize / 10       

    allocate(gbl_ipairs(ipairs_maxsize), stat = alloc_failed)

    if (alloc_failed .ne. 0) then
      call alloc_error('get_nb_list', 'ipairs array reallocation failed!');
    end if

#ifdef PAIRLST_DBG
    write(mdout, '(/,a,5x,a,i8)') &
          '|', 'Nonbonded Pairs Reallocation:', ipairs_maxsize
#endif

  end do

#ifdef MPI
  if (block_fft .eq. 0) then
    if (is_orthog .ne. 0) then
      call claim_slab_recip_imgs(atm_cnt, fraction, pbc_box, crd_idx_tbl, &
                                 gbl_img_crd, gbl_img_qterm, gbl_img_iac, &
                                 gbl_mapped_img_cnt, gbl_mapped_img_lst, &
                                 gbl_img_atm_map, gbl_used_img_map, &
                                 gbl_used_img_cnt, gbl_used_img_lst)
    else
      call claim_slab_recip_imgs_nonorthog(atm_cnt, fraction, crd_idx_tbl, &
                                           gbl_img_crd, gbl_img_qterm, &
                                           gbl_img_iac, &
                                           gbl_mapped_img_cnt, &
                                           gbl_mapped_img_lst, &
                                           gbl_img_atm_map, &
                                           gbl_used_img_map, &
                                           gbl_used_img_cnt, gbl_used_img_lst)
    end if
  else
    if (is_orthog .ne. 0) then
      call claim_blk_recip_imgs(atm_cnt, fraction, pbc_box, crd_idx_tbl, &
                                gbl_img_crd, gbl_img_qterm, gbl_img_iac, &
                                gbl_mapped_img_cnt, gbl_mapped_img_lst, &
                                gbl_img_atm_map, gbl_used_img_map, &
                                gbl_used_img_cnt, gbl_used_img_lst)
    else
      call claim_blk_recip_imgs_nonorthog(atm_cnt, fraction, crd_idx_tbl, &
                                          gbl_img_crd, gbl_img_qterm, &
                                          gbl_img_iac, gbl_mapped_img_cnt, &    
                                          gbl_mapped_img_lst, &
                                          gbl_img_atm_map, &
                                          gbl_used_img_map, &
                                          gbl_used_img_cnt, gbl_used_img_lst)
    end if
  end if

  log_used_img_cnt = dble(gbl_used_img_cnt)

  call update_pme_time(build_list_timer)
  call update_time(nonbond_time)

  call get_send_atm_lst(atm_cnt)

  call update_time(fcve_dist_time)
  call zero_pme_time()

  if (imin .eq. 0) then
    call save_used_atom_crds(atm_cnt, crd, gbl_atm_saved_crd)
  else
    ! Minimizations require all coords for skin check...
    call save_all_atom_crds(atm_cnt, crd, gbl_atm_saved_crd)
  end if
#else
  call update_pme_time(build_list_timer)
  call update_time(nonbond_time)
#endif /* MPI */


  return

end subroutine pme_list

!*******************************************************************************
!
! Subroutine:  get_nb_energy
!
! Description:
!              
! The main routine for non bond energy (vdw and hbond) as well as direct part
! of ewald sum.  It is structured for parallelism.
!
!*******************************************************************************

subroutine get_nb_energy(img_frc, img_crd, img_qterm, eed_cub, &
                         ipairs, tranvec, need_pot_enes, need_virials, &
                         eed, evdw, ehb, eedvir, virial)

  use img_mod
  use timers_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_mod
  use parallel_dat_mod
  use prmtop_dat_mod
#ifdef DIRFRC_EFS
  use ene_frc_splines_mod
#endif /* DIRFRC_EFS */

  implicit none

! Formal arguments:

  double precision, intent(in out) :: img_frc(3, *)
  double precision, intent(in)     :: img_crd(3, *)
  double precision , intent(in)    :: img_qterm(*)
  double precision, intent(in)     :: eed_cub(*)
#ifdef DIRFRC_NOVEC
  integer                          :: ipairs(*)
#else
  integer, intent(in)              :: ipairs(*)
#endif
  double precision, intent(in)     :: tranvec(1:3, 0:17)
  logical, intent(in)              :: need_pot_enes
  logical, intent(in)              :: need_virials
  double precision, intent(out)    :: eed
  double precision, intent(out)    :: evdw
  double precision, intent(out)    :: ehb
  double precision, intent(out)    :: eedvir
  double precision, intent(out)    :: virial(3, 3)

! Local variables and parameters:

  double precision      del
  double precision      dxdr
  double precision      eedtbdns_stk
  double precision      eedvir_stk, eed_stk, evdw_stk, ehb_stk
  double precision      max_nb_cut2, es_cut2, es_cut
  double precision      x_i, y_i, z_i
  double precision      x_tran(1:3, 0:17)
  double precision      vxx, vxy, vxz, vyy, vyz, vzz
  integer               i
  integer               ipairs_idx
  integer               ntypes_stk
  integer               img_i
  integer               ee_eval_cnt
  integer               full_eval_cnt
#ifdef DIRFRC_COMTRANS
  integer               common_tran    ! flag - 1 if translation not needed
#endif /* DIRFRC_COMTRANS */
  logical               cutoffs_equal

#ifdef MPI
  if (my_img_lo .gt. my_img_hi) return
#endif /* MPI */

  ntypes_stk = ntypes

  eedvir_stk = 0.d0
  eed_stk = 0.d0
  evdw_stk = 0.d0
  ehb_stk = 0.d0

  dxdr = ew_coeff
  eedtbdns_stk = eedtbdns
  del = 1.d0 / eedtbdns_stk
  max_nb_cut2 = vdw_cutoff * vdw_cutoff
  es_cut = es_cutoff
  es_cut2 = es_cut * es_cut
  cutoffs_equal = (vdw_cutoff .eq. es_cutoff)

  vxx = 0.d0
  vxy = 0.d0
  vxz = 0.d0
  vyy = 0.d0
  vyz = 0.d0
  vzz = 0.d0

  ipairs_idx = 1

  if (need_pot_enes) then

    do img_i = my_img_lo, my_img_hi

#ifdef DIRFRC_COMTRANS
      ! Common translation (ie. no translation) flag is packed at
      ! the front of each sublist followed by the count(s) of sublist
      ! image pair entries.

      common_tran = ipairs(ipairs_idx)
      ipairs_idx = ipairs_idx + 1
#endif /* DIRFRC_COMTRANS */
    
      ! Electrostatic evaluation-only count followed by
      ! full evaluation count packed at the front of each pair sublist.

      ee_eval_cnt = ipairs(ipairs_idx)
      full_eval_cnt = ipairs(ipairs_idx + 1)
      ipairs_idx = ipairs_idx + 2

      if (ee_eval_cnt + full_eval_cnt .gt. 0) then

        x_i = img_crd(1, img_i)
        y_i = img_crd(2, img_i)
        z_i = img_crd(3, img_i)

#ifdef DIRFRC_COMTRANS
        if (common_tran .eq. 0) then
#endif /* DIRFRC_COMTRANS */
          ! We need all the translation vectors:
          do i = 0, 17
            x_tran(1, i) = tranvec(1, i) - x_i
            x_tran(2, i) = tranvec(2, i) - y_i
            x_tran(3, i) = tranvec(3, i) - z_i
          end do
#ifdef DIRFRC_COMTRANS
        else
          ! Just put the x,y,z values in the middle cell
          x_tran(1, 13) = - x_i
          x_tran(2, 13) = - y_i
          x_tran(3, 13) = - z_i
        end if
#endif /* DIRFRC_COMTRANS */

        ! We always need virials from this routine if we need energies,
        ! because the virials are used in estimating the pme error.

        if (cutoffs_equal) then
#ifdef DIRFRC_EFS
          call pairs_calc_efv(img_frc, img_crd, img_qterm, efs_tbl, &
                              eed_cub, typ_ico, ipairs(ipairs_idx), &
                              gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
#else
          call pairs_calc_efv(img_frc, img_crd, img_qterm, eed_cub, &
                              typ_ico, ipairs(ipairs_idx), &
                              gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
#endif 
        else
#ifdef DIRFRC_EFS
          call pairs_calc_efv_2cut(img_frc, img_crd, img_qterm, efs_tbl, &
                                   eed_cub, typ_ico, ipairs(ipairs_idx), &
                                   gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
#else
          call pairs_calc_efv_2cut(img_frc, img_crd, img_qterm, eed_cub, &
                                   typ_ico, ipairs(ipairs_idx), &
                                   gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
#endif 
        end if

        ipairs_idx = ipairs_idx + ee_eval_cnt + full_eval_cnt

      end if
    end do

  else ! do not need energies...

    do img_i = my_img_lo, my_img_hi

#ifdef DIRFRC_COMTRANS
      ! Common translation (ie. no translation) flag is packed at
      ! the front of each sublist followed by the count(s) of sublist
      ! image pair entries.

      common_tran = ipairs(ipairs_idx)
      ipairs_idx = ipairs_idx + 1
#endif /* DIRFRC_COMTRANS */
    
      ! Electrostatic evaluation-only count followed by
      ! full evaluation count packed at the front of each pair sublist.

      ee_eval_cnt = ipairs(ipairs_idx)
      full_eval_cnt = ipairs(ipairs_idx + 1)
      ipairs_idx = ipairs_idx + 2

      if (ee_eval_cnt + full_eval_cnt .gt. 0) then

        x_i = img_crd(1, img_i)
        y_i = img_crd(2, img_i)
        z_i = img_crd(3, img_i)

#ifdef DIRFRC_COMTRANS
        if (common_tran .eq. 0) then
#endif /* DIRFRC_COMTRANS */
          ! We need all the translation vectors:
          do i = 0, 17
            x_tran(1, i) = tranvec(1, i) - x_i
            x_tran(2, i) = tranvec(2, i) - y_i
            x_tran(3, i) = tranvec(3, i) - z_i
          end do
#ifdef DIRFRC_COMTRANS
        else
          ! Just put the x,y,z values in the middle cell
          x_tran(1, 13) = - x_i
          x_tran(2, 13) = - y_i
          x_tran(3, 13) = - z_i
        end if
#endif /* DIRFRC_COMTRANS */

        if (need_virials) then

          if (cutoffs_equal) then
#ifdef DIRFRC_EFS
            call pairs_calc_fv(img_frc, img_crd, img_qterm, fs_tbl, &
                               eed_cub, typ_ico, ipairs(ipairs_idx), &
                               gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
#else
            call pairs_calc_fv(img_frc, img_crd, img_qterm, eed_cub, &
                               typ_ico, ipairs(ipairs_idx), &
                               gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
#endif 
          else
#ifdef DIRFRC_EFS
            call pairs_calc_fv_2cut(img_frc, img_crd, img_qterm, fs_tbl, &
                                    eed_cub, typ_ico, ipairs(ipairs_idx), &
                                    gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
#else
            call pairs_calc_fv_2cut(img_frc, img_crd, img_qterm, eed_cub, &
                                    typ_ico, ipairs(ipairs_idx), &
                                    gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
#endif 
          end if

        else

          if (cutoffs_equal) then
#ifdef DIRFRC_EFS
            call pairs_calc_f(img_frc, img_crd, img_qterm, fs_tbl, &
                              eed_cub, typ_ico, ipairs(ipairs_idx), &
                              gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
#else
            call pairs_calc_f(img_frc, img_crd, img_qterm, eed_cub, &
                              typ_ico, ipairs(ipairs_idx), &
                              gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
#endif 
          else
#ifdef DIRFRC_EFS
            call pairs_calc_f_2cut(img_frc, img_crd, img_qterm, fs_tbl, &
                                   eed_cub, typ_ico, ipairs(ipairs_idx), &
                                   gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
#else
            call pairs_calc_f_2cut(img_frc, img_crd, img_qterm, eed_cub, &
                                   typ_ico, ipairs(ipairs_idx), &
                                   gbl_img_iac, gbl_cn1, gbl_cn2, x_tran)
#endif 
          end if

        end if

        ipairs_idx = ipairs_idx + ee_eval_cnt + full_eval_cnt

      end if
    end do

  end if

  ! Save the energies:
                                                                                
  eedvir = eedvir_stk
  eed = eed_stk
  evdw = evdw_stk
  ehb = ehb_stk

  ! Save the virials.

  virial(1, 1) = vxx
  virial(1, 2) = vxy
  virial(2, 1) = vxy
  virial(1, 3) = vxz
  virial(3, 1) = vxz
  virial(2, 2) = vyy
  virial(2, 3) = vyz
  virial(3, 2) = vyz
  virial(3, 3) = vzz

  return

contains

#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#include "pairs_calc_vec.i"
#include "pairs_calc_novec.i"
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#include "pairs_calc_vec.i"
#include "pairs_calc_novec.i"
#undef NEED_VIR
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_F
#include "pairs_calc_vec.i"
#include "pairs_calc_novec.i"
#undef BUILD_PAIRS_CALC_F

#define BUILD_PAIRS_CALC_EFV
#define NEED_ENE
#define NEED_VIR
#define BUILD_PAIRS_CALC_VEC_2CUT
#include "pairs_calc_vec.i"
#undef BUILD_PAIRS_CALC_VEC_2CUT
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include "pairs_calc_novec.i"
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef NEED_VIR
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_EFV

#define BUILD_PAIRS_CALC_FV
#define NEED_VIR
#define BUILD_PAIRS_CALC_VEC_2CUT
#include "pairs_calc_vec.i"
#undef BUILD_PAIRS_CALC_VEC_2CUT
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include "pairs_calc_novec.i"
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef NEED_ENE
#undef BUILD_PAIRS_CALC_FV

#define BUILD_PAIRS_CALC_F
#define BUILD_PAIRS_CALC_VEC_2CUT
#include "pairs_calc_vec.i"
#undef BUILD_PAIRS_CALC_VEC_2CUT
#define BUILD_PAIRS_CALC_NOVEC_2CUT
#include "pairs_calc_novec.i"
#undef BUILD_PAIRS_CALC_NOVEC_2CUT
#undef BUILD_PAIRS_CALC_F

end subroutine get_nb_energy

end module pme_direct_mod
