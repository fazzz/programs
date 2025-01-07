#include "copyright.i"

!*******************************************************************************
!
! Module:  nb_pairlist_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module nb_pairlist_mod

  implicit none

! Data that is not broadcast:

  double precision, save, allocatable   :: gbl_atm_saved_crd(:,:)
  double precision, save, allocatable   :: gbl_saved_imgcrd(:,:)
  integer, allocatable, save            :: gbl_ipairs(:)
  double precision, save                :: gbl_saved_box(3)
  integer, save                         :: ipairs_maxsize

contains

!*******************************************************************************
!
! Subroutine:  alloc_nb_pairlist_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_nb_pairlist_mem(atm_cnt, cutlist, num_ints, num_reals)

  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  double precision, intent(in)  :: cutlist

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer                       :: alloc_failed
  integer, parameter            :: ipairs_size_min = 1000
  double precision, parameter   :: ipairs_size_coef = 0.225d0 

  allocate(gbl_atm_saved_crd(3, atm_cnt), &
           gbl_saved_imgcrd(3, atm_cnt), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(gbl_atm_saved_crd) + &
                          size(gbl_saved_imgcrd)

! Allocate pairs list:

! The following ipairs_maxsize calc is a crude heuristic assuming that the
! number of nonbonded pairs roughly scales with the cutoff volume and
! the number of atoms in the system.  If MPI is running, the total is
! divided up among the processes.  The 2 or 3 is for the counters and flags at
! the head of each atom sublist.

#ifdef DIRFRC_COMTRANS
  ipairs_maxsize = int((ipairs_size_coef * cutlist**3 + 3.d0) * dble(atm_cnt))
#else
  ipairs_maxsize = int((ipairs_size_coef * cutlist**3 + 2.d0) * dble(atm_cnt))
#endif

#ifdef MPI
  if (numtasks .gt. 0) ipairs_maxsize = ipairs_maxsize / numtasks
#endif /* MPI */

  if (ipairs_maxsize .lt. ipairs_size_min) ipairs_maxsize = ipairs_size_min

  allocate(gbl_ipairs(ipairs_maxsize), stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  return

end subroutine alloc_nb_pairlist_mem
  
#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_nb_pairlist_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_nb_pairlist_dat(atm_cnt, cutlist)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  double precision, intent(in)  :: cutlist

! Local variables:

  integer               :: num_ints, num_reals  ! returned values discarded
  
  if (.not. master) then
    num_ints = 0
    num_reals = 0
    call alloc_nb_pairlist_mem(atm_cnt, cutlist, num_ints, num_reals)
  end if

  ! The allocated data is not initialized from the master node.

  return

end subroutine bcast_nb_pairlist_dat
#endif

!*******************************************************************************
!
! Subroutine:   map_pairlist_imgs
!
! Description:  TBS
!
!*******************************************************************************

#ifdef MPI
subroutine map_pairlist_imgs(atm_cnt, flat_cit, fraction, charge, iac, &
                             img_crd, img_qterm, img_iac, mapped_img_cnt, &
                             mapped_img_lst, img_atm_map)
#else
subroutine map_pairlist_imgs(atm_cnt, fraction, charge, iac, &
                             img_crd, img_qterm, img_iac, img_atm_map)
#endif

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
#ifdef MPI
  type(cit_tbl_rec)     :: flat_cit(0 : cit_tbl_x_dim * &
                                        cit_tbl_y_dim * &
                                        cit_tbl_z_dim - 1)
#endif /* MPI */
  double precision      :: fraction(3, atm_cnt)
  double precision      :: charge(atm_cnt)
  integer               :: iac(atm_cnt)
#ifdef MPI
  integer               :: mapped_img_cnt
  integer               :: mapped_img_lst(*)
#endif /* MPI */
  double precision      :: img_crd(3, atm_cnt)
  double precision      :: img_qterm(atm_cnt)
  integer               :: img_iac(atm_cnt)
  integer               :: img_atm_map(atm_cnt)

! Local variables:

! "bkt" refers to a bucket or cell in the flat cit table.
! "bkts" refers to the bucket mapping tables.

  double precision      :: f1, f2, f3
  double precision      :: x_box, y_box, z_box
  double precision      :: ucell_stk(3, 3)
#ifdef MPI
  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
#endif /* MPI */
  integer               :: atm_idx, img_idx
  integer               :: is_orthog_stk
#ifdef MPI
  integer               :: lst_idx
  integer               :: i_bkt_x_idx, i_bkt_y_idx, i_bkt_z_idx
  integer               :: x_bkts_lo, x_bkts_hi
  integer               :: y_bkts_lo, y_bkts_hi
  integer               :: z_bkts_lo, z_bkts_hi
  integer               :: img_i_lo, atm_i_lo
  integer               :: i_bkt_lo, i_bkt_hi
  integer               :: i_bkt, cur_bkt, yz_bkt, z_bkt
  integer               :: x_bkts_idx, y_bkts_idx, z_bkts_idx
  integer               :: x_bkts(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: y_bkts(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: z_bkts(0 : cit_tbl_z_dim * 2 - 1)
#endif /* MPI */

  x_box = pbc_box(1)
  y_box = pbc_box(2)
  z_box = pbc_box(3)

  ucell_stk(:,:) = ucell(:,:)

  is_orthog_stk = is_orthog

#ifdef MPI

  do lst_idx = 1, mapped_img_cnt
    img_idx = mapped_img_lst(lst_idx)
    img_iac(img_idx) = 0
  end do
  mapped_img_cnt = 0

  if (my_img_lo .gt. my_img_hi) return

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

! Set up the bucket mapping arrays.  Redoing this is only necessary if the
! box dimensions change, but it is cheap and there is a locality benefit
! to having it on the stack.

  call setup_cit_tbl_bkts(x_bkts, y_bkts, z_bkts)

! Get the low and high flat cit bucket indexes for this task:

  call get_flat_cit_idx(my_img_lo, img_atm_map, fraction, i_bkt_lo)
  call get_flat_cit_idx(my_img_hi, img_atm_map, fraction, i_bkt_hi)

! Map all the images you will reference in building the pairlist.
! Note that "mapping" is different than claiming the images as "used".
! When mapped, this just means all the image's data structures are valid.

  do i_bkt = i_bkt_lo, i_bkt_hi

    img_i_lo = flat_cit(i_bkt)%img_lo
    if (img_i_lo .eq. 0) cycle

    ! Get the x,y,z 3d cit indexes to be able to get the bucket ranges
    ! to search:

    atm_i_lo = img_atm_map(img_i_lo)
    i_bkt_x_idx = int(fraction(1, atm_i_lo) * scale_fac_x)
    i_bkt_y_idx = int(fraction(2, atm_i_lo) * scale_fac_y)
    i_bkt_z_idx = int(fraction(3, atm_i_lo) * scale_fac_z)

    ! Determine the bucket ranges that need to be searched:

    x_bkts_lo = i_bkt_x_idx + cit_tbl_x_dim - cit_bkt_delta
    x_bkts_hi = i_bkt_x_idx + cit_tbl_x_dim + cit_bkt_delta
    y_bkts_lo = i_bkt_y_idx + cit_tbl_y_dim - cit_bkt_delta
    y_bkts_hi = i_bkt_y_idx + cit_tbl_y_dim + cit_bkt_delta
    z_bkts_lo = i_bkt_z_idx + cit_tbl_z_dim - cit_bkt_delta
    z_bkts_hi = i_bkt_z_idx + cit_tbl_z_dim
                                        
    z_loop: &
    do z_bkts_idx = z_bkts_lo, z_bkts_hi
      z_bkt = z_bkts(z_bkts_idx)
      y_loop: &
      do y_bkts_idx = y_bkts_lo, y_bkts_hi
        yz_bkt = z_bkt + y_bkts(y_bkts_idx)
        x_loop: &
        do x_bkts_idx = x_bkts_lo, x_bkts_hi
          cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
          img_i_lo = flat_cit(cur_bkt)%img_lo
          if (img_i_lo .eq. 0) cycle
          if (img_iac(img_i_lo) .ne. 0) cycle
          if (is_orthog_stk .ne. 0) then
            do img_idx = img_i_lo, flat_cit(cur_bkt)%img_hi
              atm_idx = img_atm_map(img_idx)
              img_crd(1, img_idx) = fraction(1, atm_idx) * x_box
              img_crd(2, img_idx) = fraction(2, atm_idx) * y_box
              img_crd(3, img_idx) = fraction(3, atm_idx) * z_box
              img_qterm(img_idx) = charge(atm_idx)
              img_iac(img_idx) = iac(atm_idx)
              mapped_img_cnt = mapped_img_cnt + 1
              mapped_img_lst(mapped_img_cnt) = img_idx
            end do
          else
            do img_idx = img_i_lo, flat_cit(cur_bkt)%img_hi
              atm_idx = img_atm_map(img_idx)
              f1 = fraction(1, atm_idx)
              f2 = fraction(2, atm_idx)
              f3 = fraction(3, atm_idx)
              ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
              ! we can simplify the expression in this critical inner loop
              img_crd(1, img_idx) = f1 * ucell_stk(1, 1) + &
                                    f2 * ucell_stk(1, 2) + &
                                    f3 * ucell_stk(1, 3)
              img_crd(2, img_idx) = f2 * ucell_stk(2, 2) + &
                                    f3 * ucell_stk(2, 3)
              img_crd(3, img_idx) = f3 * ucell_stk(3, 3)
              img_qterm(img_idx) = charge(atm_idx)
              img_iac(img_idx) = iac(atm_idx)
              mapped_img_cnt = mapped_img_cnt + 1
              mapped_img_lst(mapped_img_cnt) = img_idx
            end do
          end if
          if (cur_bkt .eq. i_bkt) exit z_loop ! Done with cit buckets
                                             ! image i claims pairs from.
        end do x_loop
      end do y_loop
    end do z_loop

  end do

#else

! For the uniprocessor, just map it all...

  if (is_orthog_stk .ne. 0) then
    do img_idx = 1, atm_cnt
      atm_idx = img_atm_map(img_idx)
      img_crd(1, img_idx) = fraction(1, atm_idx) * x_box
      img_crd(2, img_idx) = fraction(2, atm_idx) * y_box
      img_crd(3, img_idx) = fraction(3, atm_idx) * z_box
      img_qterm(img_idx) = charge(atm_idx)
      img_iac(img_idx) = iac(atm_idx)
    end do
  else
    do img_idx = 1, atm_cnt
      atm_idx = img_atm_map(img_idx)
      f1 = fraction(1, atm_idx)
      f2 = fraction(2, atm_idx)
      f3 = fraction(3, atm_idx)
      ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
      ! we can simplify the expression in this critical inner loop
      img_crd(1, img_idx) = f1 * ucell_stk(1, 1) + &
                            f2 * ucell_stk(1, 2) + &
                            f3 * ucell_stk(1, 3)
      img_crd(2, img_idx) = f2 * ucell_stk(2, 2) + &
                            f3 * ucell_stk(2, 3)
      img_crd(3, img_idx) = f3 * ucell_stk(3, 3)
      img_qterm(img_idx) = charge(atm_idx)
      img_iac(img_idx) = iac(atm_idx)
    end do
  end if

#endif /* MPI */

  return

end subroutine map_pairlist_imgs

!*******************************************************************************
!
! Subroutine:   get_nb_list
!
! Description:  Main routine for pairlist building.  
!
!*******************************************************************************

subroutine get_nb_list(atm_cnt, flat_cit, img_crd, img_atm_map, &
#ifdef MPI
                       used_img_map, used_img_cnt, used_img_lst, &
#endif
                       ico, fraction, tranvec, atm_maskdata, atm_mask, &
                       atm_img_map, excl_img_flags, &
                       img_iac, igroup, ntypes, ibelly, &
                       es_cutoff, vdw_cutoff, skinnb, verbose, ifail)

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  type(cit_tbl_rec)     :: flat_cit(0 : cit_tbl_x_dim * &
                                        cit_tbl_y_dim * &
                                        cit_tbl_z_dim - 1)

  double precision      :: img_crd(3, atm_cnt)
  integer               :: img_atm_map(atm_cnt)
#ifdef MPI
  integer(byte)         :: used_img_map(atm_cnt)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
#endif
  integer               :: ico(*)
  double precision      :: fraction(3, atm_cnt)
  double precision      :: tranvec(1:3, 0:17)
  type(listdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)
  integer               :: atm_img_map(atm_cnt)
  integer               :: excl_img_flags(atm_cnt)
  integer               :: img_iac(atm_cnt)
  integer               :: igroup(atm_cnt)
  integer               :: ntypes
  integer               :: ibelly
  double precision      :: es_cutoff
  double precision      :: vdw_cutoff
  double precision      :: skinnb
  integer               :: verbose
  integer               :: ifail

! Local variables:

  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
  double precision      :: cutlist_sq
  double precision      :: x_i, y_i, z_i
  integer               :: atm_i, img_i
  integer               :: i
  integer               :: atm_mask_idx
  integer               :: num_packed
  integer               :: i_bkt_x_idx, i_bkt_y_idx, i_bkt_z_idx
  integer               :: x_bkts_lo, x_bkts_hi
  integer               :: y_bkts_lo, y_bkts_hi
  integer               :: z_bkts_lo, z_bkts_hi
  integer               :: i_bkt, saved_i_bkt_img_hi
  integer               :: img_i_lo, img_i_hi
  integer               :: i_bkt_lo, i_bkt_hi
#ifdef DIRFRC_COMTRANS
  logical               :: common_tran
#endif /* DIRFRC_COMTRANS */
  logical               :: dont_skip_belly_pairs
  integer               :: iaci, ntypes_stk
  integer               :: ee_eval_cnt
  integer               :: full_eval_cnt
  integer               :: is_orthog_stk

! Variables use in the included files.  We would put this stuff in 
! contained subroutines, but it really whacks performance for itanium/ifort.

  double precision      :: dx, dy, dz
  integer               :: cur_bkt, yz_bkt, z_bkt
  integer               :: cur_tran, yz_tran, z_tran
  integer               :: img_j
  integer               :: x_bkts_idx, y_bkts_idx, z_bkts_idx
  double precision      :: r_sq
  integer               :: enc_img
  double precision      :: i_tranvec(1:3, 0:17)

! Variables used for double cutoffs:

  double precision      :: es_cut_sq
  logical               :: cutoffs_equal

  integer               :: x_bkts(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: x_trans(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: y_bkts(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: y_trans(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: z_bkts(0 : cit_tbl_z_dim * 2 - 1)
  integer               :: z_trans(0 : cit_tbl_z_dim * 2 - 1)

  integer               :: total_pairs  ! DBG

  integer               :: img_j_ee_eval(atm_cnt)
  integer               :: img_j_full_eval(atm_cnt)

! "bkt" refers to a bucket or cell in the flat cit table.
! "bkts" refers to the bucket mapping tables.
! "trans" refers to the translation mapping tables.

#ifdef MPI
  if (my_img_lo .gt. my_img_hi) return
#endif /* MPI */
  num_packed = 0        ! For "regular" pairlist
  total_pairs = 0       ! DBG

  dont_skip_belly_pairs = ibelly .eq. 0
  ntypes_stk = ntypes
  is_orthog_stk = is_orthog

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

  es_cut_sq = (es_cutoff + skinnb) * (es_cutoff + skinnb)
  cutoffs_equal = (es_cutoff .eq. vdw_cutoff)
  cutlist_sq = (vdw_cutoff + skinnb) * (vdw_cutoff + skinnb)

  ! Set up the bucket and translation index mapping arrays.  We really only need
  ! to do this each time if the cit table dimensions are changing, but there
  ! may be advantages to having only local variables anyway.

  call setup_cit_tbl_bkts(x_bkts, y_bkts, z_bkts, x_trans, y_trans, z_trans)

  ! Get the low and high flat cit bucket indexes for this task:

#ifdef MPI
  call get_flat_cit_idx(my_img_lo, img_atm_map, fraction, i_bkt_lo)
  call get_flat_cit_idx(my_img_hi, img_atm_map, fraction, i_bkt_hi)
#else
  i_bkt_lo = 0
  i_bkt_hi = cit_tbl_x_dim * cit_tbl_y_dim * cit_tbl_z_dim - 1
#endif /* MPI */

  ! Loop over the atoms you own, finding their nonbonded partners...
  ! img_i is the atom for which you are finding the partner atoms, img_j

  do i_bkt = i_bkt_lo, i_bkt_hi

    img_i_lo = flat_cit(i_bkt)%img_lo
    img_i_hi = flat_cit(i_bkt)%img_hi

    if (img_i_lo .eq. 0) cycle
    saved_i_bkt_img_hi = img_i_hi
#ifdef MPI
    if (img_i_lo .lt. my_img_lo) img_i_lo = my_img_lo
    if (img_i_hi .gt. my_img_hi) img_i_hi = my_img_hi
#endif /* MPI */

    ! Get the x,y,z 3d cit indexes to be able to get the bucket ranges
    ! to search:

    atm_i = img_atm_map(img_i_lo)
    i_bkt_x_idx = int(fraction(1, atm_i) * scale_fac_x)
    i_bkt_y_idx = int(fraction(2, atm_i) * scale_fac_y)
    i_bkt_z_idx = int(fraction(3, atm_i) * scale_fac_z)

    ! Determine the bucket ranges that need to be searched:

    x_bkts_lo = i_bkt_x_idx + cit_tbl_x_dim - cit_bkt_delta
    x_bkts_hi = i_bkt_x_idx + cit_tbl_x_dim + cit_bkt_delta
    y_bkts_lo = i_bkt_y_idx + cit_tbl_y_dim - cit_bkt_delta
    y_bkts_hi = i_bkt_y_idx + cit_tbl_y_dim + cit_bkt_delta
    z_bkts_lo = i_bkt_z_idx + cit_tbl_z_dim - cit_bkt_delta
    z_bkts_hi = i_bkt_z_idx + cit_tbl_z_dim
                                        
#ifdef DIRFRC_COMTRANS
    common_tran = x_trans(x_bkts_lo) .eq. 1 .and. &
                  x_trans(x_bkts_hi) .eq. 1 .and. &
                  y_trans(y_bkts_lo) .eq. 3 .and. &
                  y_trans(y_bkts_hi) .eq. 3 .and. &
                  z_trans(z_bkts_lo) .eq. 9 .and. &
                  z_trans(z_bkts_hi) .eq. 9
#endif /* DIRFRC_COMTRANS */

    do img_i = img_i_lo, img_i_hi

      flat_cit(i_bkt)%img_hi = img_i - 1

      ee_eval_cnt = 0
      full_eval_cnt = 0
      iaci = ntypes_stk * (img_iac(img_i) - 1)

      atm_i = img_atm_map(img_i)

      ! Mark excluded j images:

      atm_mask_idx = atm_maskdata(atm_i)%offset
      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 1
      end do

      ! These are imaged coordinates:

      x_i = img_crd(1, img_i)
      y_i = img_crd(2, img_i)
      z_i = img_crd(3, img_i)

#ifdef DIRFRC_COMTRANS

      if (cutoffs_equal) then

        if (common_tran) then

          z_loop1: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop1: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop1: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = img_j
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop1 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop1
            end do y_loop1
          end do z_loop1

        else    ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i 
            i_tranvec(2, i) = tranvec(2, i) - y_i 
            i_tranvec(3, i) = tranvec(3, i) - z_i 
          end do

          z_loop2: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop2: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop2: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      enc_img = ior(img_j, ishft(cur_tran, 27))
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = enc_img
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = enc_img
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop2 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop2
            end do y_loop2
          end do z_loop2

        end if

      else  ! cutoffs not equal

        if (common_tran) then

          z_loop3: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop3: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop3: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = img_j
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop3 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop3
            end do y_loop3
          end do z_loop3

        else  ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i 
            i_tranvec(2, i) = tranvec(2, i) - y_i 
            i_tranvec(3, i) = tranvec(3, i) - z_i 
          end do

          z_loop4: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop4: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop4: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = &
                            ior(img_j,ishft(cur_tran, 27))
#ifdef MPI
                          if (used_img_map(img_j) .eq. 0) then
                            used_img_map(img_j) = 1
                            used_img_cnt = used_img_cnt + 1
                            used_img_lst(used_img_cnt) = img_j
                          end if
#endif /* MPI */
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = &
                          ior(img_j,ishft(cur_tran,27))
#ifdef MPI
                        if (used_img_map(img_j) .eq. 0) then
                          used_img_map(img_j) = 1
                          used_img_cnt = used_img_cnt + 1
                          used_img_lst(used_img_cnt) = img_j
                        end if
#endif /* MPI */
                      end if
                    end if
                  end if
                end do

                if (cur_bkt .eq. i_bkt) exit z_loop4 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop4
            end do y_loop4
          end do z_loop4
        end if
      end if
#else  /* NOT DIRFRC_COMTRANS */
      if (cutoffs_equal) then

        do i = 0, 17
          i_tranvec(1, i) = tranvec(1, i) - x_i 
          i_tranvec(2, i) = tranvec(2, i) - y_i 
          i_tranvec(3, i) = tranvec(3, i) - z_i 
        end do

        z_loop2: &
        do z_bkts_idx = z_bkts_lo, z_bkts_hi
          z_bkt = z_bkts(z_bkts_idx)
          z_tran = z_trans(z_bkts_idx)
          y_loop2: &
          do y_bkts_idx = y_bkts_lo, y_bkts_hi
            yz_bkt = z_bkt + y_bkts(y_bkts_idx)
            yz_tran = z_tran + y_trans(y_bkts_idx)
            x_loop2: &
            do x_bkts_idx = x_bkts_lo, x_bkts_hi
              cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
              cur_tran = x_trans(x_bkts_idx) + yz_tran
              do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                  if (excl_img_flags(img_j) .eq. 0) then
                    enc_img = ior(img_j, ishft(cur_tran, 27))
                    if (ico(iaci + img_iac(img_j)) .eq. 0) then
                      ee_eval_cnt = ee_eval_cnt + 1
                      img_j_ee_eval(ee_eval_cnt) = enc_img
                    else
                      full_eval_cnt = full_eval_cnt + 1
                      img_j_full_eval(full_eval_cnt) = enc_img
                    end if
#ifdef MPI
                    if (used_img_map(img_j) .eq. 0) then
                      used_img_map(img_j) = 1
                      used_img_cnt = used_img_cnt + 1
                      used_img_lst(used_img_cnt) = img_j
                    end if
#endif /* MPI */
                  end if
                end if
              end do
              if (cur_bkt .eq. i_bkt) exit z_loop2 ! Done with cit buckets
                                                   ! image i claims pairs from.
            end do x_loop2
          end do y_loop2
        end do z_loop2

      else      ! cutoffs not equal

        do i = 0, 17
          i_tranvec(1, i) = tranvec(1, i) - x_i 
          i_tranvec(2, i) = tranvec(2, i) - y_i 
          i_tranvec(3, i) = tranvec(3, i) - z_i 
        end do

        z_loop4: &
        do z_bkts_idx = z_bkts_lo, z_bkts_hi
          z_bkt = z_bkts(z_bkts_idx)
          z_tran = z_trans(z_bkts_idx)
          y_loop4: &
          do y_bkts_idx = y_bkts_lo, y_bkts_hi
            yz_bkt = z_bkt + y_bkts(y_bkts_idx)
            yz_tran = z_tran + y_trans(y_bkts_idx)
            x_loop4: &
            do x_bkts_idx = x_bkts_lo, x_bkts_hi
              cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
              cur_tran = x_trans(x_bkts_idx) + yz_tran
              do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                r_sq = dx * dx + dy * dy + dz * dz
                if (r_sq .lt. cutlist_sq) then
                  if (excl_img_flags(img_j) .eq. 0) then
                    if (ico(iaci + img_iac(img_j)) .eq. 0) then
                      if (r_sq .lt. es_cut_sq) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = &
                          ior(img_j,ishft(cur_tran, 27))
#ifdef MPI
                        if (used_img_map(img_j) .eq. 0) then
                          used_img_map(img_j) = 1
                          used_img_cnt = used_img_cnt + 1
                          used_img_lst(used_img_cnt) = img_j
                        end if
#endif /* MPI */
                      end if
                    else
                      full_eval_cnt = full_eval_cnt + 1
                      img_j_full_eval(full_eval_cnt) = &
                        ior(img_j,ishft(cur_tran,27))
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end if
              end do
              if (cur_bkt .eq. i_bkt) exit z_loop4 ! Done with cit buckets
                                                   ! image i claims pairs from.
            end do x_loop4
          end do y_loop4
        end do z_loop4

      end if
#endif /* DIRFRC_COMTRANS */

      total_pairs = total_pairs + full_eval_cnt + ee_eval_cnt       ! DBG

      if (dont_skip_belly_pairs) then
        call pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                          full_eval_cnt, img_j_full_eval, &
                          gbl_ipairs, num_packed)
      else
        call pack_nb_list_skip_belly_pairs(ee_eval_cnt, img_j_ee_eval, &
                                           full_eval_cnt, img_j_full_eval, &
                                           gbl_ipairs, igroup, img_atm_map, &
                                           num_packed)
      end if

      ! Clear excluded j images flags:

      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 0
      end do

      if (ifail .eq. 1) then
        flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi
        return
      end if

    end do ! end of images in i_bkt

    flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi

  end do ! end of all cit buckets with i images owned by this task

! write(0, *) 'DBG: total_pairs, listtot = ', &
!             total_pairs, num_packed
! write(0, *) 'DBG: avg packing = ', &
!             dble(total_pairs)/dble(num_packed)

  if (verbose .gt. 0) then
    if (master) write(mdout, *) 'listtot = ', num_packed
  end if

  return

contains

#include "pack_nb_list_dflt.i"

end subroutine get_nb_list

!*******************************************************************************
!
! Subroutine:   get_nb_ips_list
!
! Description:  Main routine for pairlist building.  
!
!*******************************************************************************

subroutine get_nb_ips_list(atm_cnt, flat_cit, img_crd, img_atm_map, &
#ifdef MPI
                       used_img_map, used_img_cnt, used_img_lst, &
                       used_img_ips_map, used_img_ips_cnt, used_img_ips_lst, &
#endif
                       ico, fraction, tranvec, atm_maskdata, atm_mask, &
                       atm_img_map, excl_img_flags, &
                       img_iac, igroup, ntypes, ibelly, &
                       es_cutoff, vdw_cutoff, skinnb, verbose, ifail)

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  type(cit_tbl_rec)     :: flat_cit(0 : cit_tbl_x_dim * &
                                        cit_tbl_y_dim * &
                                        cit_tbl_z_dim - 1)

  double precision      :: img_crd(3, atm_cnt)
  integer               :: img_atm_map(atm_cnt)
#ifdef MPI
  integer(byte)         :: used_img_map(atm_cnt)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
  integer(byte)         :: used_img_ips_map(atm_cnt)
  integer               :: used_img_ips_cnt
  integer               :: used_img_ips_lst(*)
#endif
  integer               :: ico(*)
  double precision      :: fraction(3, atm_cnt)
  double precision      :: tranvec(1:3, 0:17)
  type(listdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)
  integer               :: atm_img_map(atm_cnt)
  integer               :: excl_img_flags(atm_cnt)
  integer               :: img_iac(atm_cnt)
  integer               :: igroup(atm_cnt)
  integer               :: ntypes
  integer               :: ibelly
  double precision      :: es_cutoff
  double precision      :: vdw_cutoff
  double precision      :: skinnb
  integer               :: verbose
  integer               :: ifail

! Local variables:

  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
  double precision      :: cutlist_sq
  double precision      :: x_i, y_i, z_i
  integer               :: atm_i, img_i
  integer               :: i
  integer               :: atm_mask_idx
  integer               :: num_packed
  integer               :: i_bkt_x_idx, i_bkt_y_idx, i_bkt_z_idx
  integer               :: x_bkts_lo, x_bkts_hi
  integer               :: y_bkts_lo, y_bkts_hi
  integer               :: z_bkts_lo, z_bkts_hi
  integer               :: i_bkt, saved_i_bkt_img_hi
  integer               :: img_i_lo, img_i_hi
  integer               :: i_bkt_lo, i_bkt_hi
#ifdef DIRFRC_COMTRANS
  logical               :: common_tran
#endif /* DIRFRC_COMTRANS */
  logical               :: dont_skip_belly_pairs
  integer               :: iaci, ntypes_stk
  integer               :: ee_eval_cnt
  integer               :: full_eval_cnt
  integer               :: is_orthog_stk

! Variables use in the included files.  We would put this stuff in 
! contained subroutines, but it really whacks performance for itanium/ifort.

  double precision      :: dx, dy, dz
  integer               :: cur_bkt, yz_bkt, z_bkt
  integer               :: cur_tran, yz_tran, z_tran
  integer               :: img_j
  integer               :: x_bkts_idx, y_bkts_idx, z_bkts_idx
  double precision      :: r_sq
  integer               :: enc_img
  double precision      :: i_tranvec(1:3, 0:17)

! Variables used for double cutoffs:

  double precision      :: es_cut_sq
  logical               :: cutoffs_equal

  integer               :: x_bkts(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: x_trans(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: y_bkts(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: y_trans(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: z_bkts(0 : cit_tbl_z_dim * 2 - 1)
  integer               :: z_trans(0 : cit_tbl_z_dim * 2 - 1)

  integer               :: total_pairs  ! DBG

  integer               :: img_j_ee_eval(atm_cnt)
  integer               :: img_j_full_eval(atm_cnt)

! "bkt" refers to a bucket or cell in the flat cit table.
! "bkts" refers to the bucket mapping tables.
! "trans" refers to the translation mapping tables.

#ifdef MPI
  if (my_img_lo .gt. my_img_hi) return
#endif /* MPI */
  num_packed = 0        ! For "regular" pairlist
  total_pairs = 0       ! DBG

  dont_skip_belly_pairs = ibelly .eq. 0
  ntypes_stk = ntypes
  is_orthog_stk = is_orthog

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

  es_cut_sq = (es_cutoff + skinnb) * (es_cutoff + skinnb)
  cutoffs_equal = (es_cutoff .eq. vdw_cutoff)
  cutlist_sq = (vdw_cutoff + skinnb) * (vdw_cutoff + skinnb)

  ! Set up the bucket and translation index mapping arrays.  We really only need
  ! to do this each time if the cit table dimensions are changing, but there
  ! may be advantages to having only local variables anyway.

  call setup_cit_tbl_bkts(x_bkts, y_bkts, z_bkts, x_trans, y_trans, z_trans)

  ! Get the low and high flat cit bucket indexes for this task:

#ifdef MPI
  call get_flat_cit_idx(my_img_lo, img_atm_map, fraction, i_bkt_lo)
  call get_flat_cit_idx(my_img_hi, img_atm_map, fraction, i_bkt_hi)
#else
  i_bkt_lo = 0
  i_bkt_hi = cit_tbl_x_dim * cit_tbl_y_dim * cit_tbl_z_dim - 1
#endif /* MPI */

  ! Loop over the atoms you own, finding their nonbonded partners...
  ! img_i is the atom for which you are finding the partner atoms, img_j

  do i_bkt = i_bkt_lo, i_bkt_hi

    img_i_lo = flat_cit(i_bkt)%img_lo
    img_i_hi = flat_cit(i_bkt)%img_hi

    if (img_i_lo .eq. 0) cycle
    saved_i_bkt_img_hi = img_i_hi
#ifdef MPI
    if (img_i_lo .lt. my_img_lo) img_i_lo = my_img_lo
    if (img_i_hi .gt. my_img_hi) img_i_hi = my_img_hi
#endif /* MPI */

    ! Get the x,y,z 3d cit indexes to be able to get the bucket ranges
    ! to search:

    atm_i = img_atm_map(img_i_lo)
    i_bkt_x_idx = int(fraction(1, atm_i) * scale_fac_x)
    i_bkt_y_idx = int(fraction(2, atm_i) * scale_fac_y)
    i_bkt_z_idx = int(fraction(3, atm_i) * scale_fac_z)

    ! Determine the bucket ranges that need to be searched:

    x_bkts_lo = i_bkt_x_idx + cit_tbl_x_dim - cit_bkt_delta
    x_bkts_hi = i_bkt_x_idx + cit_tbl_x_dim + cit_bkt_delta
    y_bkts_lo = i_bkt_y_idx + cit_tbl_y_dim - cit_bkt_delta
    y_bkts_hi = i_bkt_y_idx + cit_tbl_y_dim + cit_bkt_delta
    z_bkts_lo = i_bkt_z_idx + cit_tbl_z_dim - cit_bkt_delta
    z_bkts_hi = i_bkt_z_idx + cit_tbl_z_dim
                                        
#ifdef DIRFRC_COMTRANS
    common_tran = x_trans(x_bkts_lo) .eq. 1 .and. &
                  x_trans(x_bkts_hi) .eq. 1 .and. &
                  y_trans(y_bkts_lo) .eq. 3 .and. &
                  y_trans(y_bkts_hi) .eq. 3 .and. &
                  z_trans(z_bkts_lo) .eq. 9 .and. &
                  z_trans(z_bkts_hi) .eq. 9
#endif /* DIRFRC_COMTRANS */

    do img_i = img_i_lo, img_i_hi

      flat_cit(i_bkt)%img_hi = img_i - 1

      ee_eval_cnt = 0
      full_eval_cnt = 0
      iaci = ntypes_stk * (img_iac(img_i) - 1)

      atm_i = img_atm_map(img_i)

      ! Mark excluded j images:

      atm_mask_idx = atm_maskdata(atm_i)%offset
      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 1
      end do

      ! These are imaged coordinates:

      x_i = img_crd(1, img_i)
      y_i = img_crd(2, img_i)
      z_i = img_crd(3, img_i)

#ifdef DIRFRC_COMTRANS

      if (cutoffs_equal) then

        if (common_tran) then

          z_loop1: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop1: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop1: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = img_j
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
                    else
                      if (used_img_ips_map(img_j) .eq. 0) then
                        used_img_ips_map(img_j) = 1
                        used_img_ips_cnt = used_img_ips_cnt + 1
                        used_img_ips_lst(used_img_ips_cnt) = img_j                      
                      end if
#endif /* MPI */
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop1 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop1
            end do y_loop1
          end do z_loop1

        else    ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i 
            i_tranvec(2, i) = tranvec(2, i) - y_i 
            i_tranvec(3, i) = tranvec(3, i) - z_i 
          end do

          z_loop2: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop2: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop2: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      enc_img = ior(img_j, ishft(cur_tran, 27))
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = enc_img
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = enc_img
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
                    else 
                      if (used_img_ips_map(img_j) .eq. 0) then
                        used_img_ips_map(img_j) = 1
                        used_img_ips_cnt = used_img_ips_cnt + 1
                        used_img_ips_lst(used_img_ips_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop2 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop2
            end do y_loop2
          end do z_loop2

        end if

      else  ! cutoffs not equal

        if (common_tran) then

          z_loop3: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop3: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop3: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = img_j
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
                    else
                      if (used_img_ips_map(img_j) .eq. 0) then
                        used_img_ips_map(img_j) = 1
                        used_img_ips_cnt = used_img_ips_cnt + 1
                        used_img_ips_lst(used_img_ips_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop3 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop3
            end do y_loop3
          end do z_loop3

        else  ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i 
            i_tranvec(2, i) = tranvec(2, i) - y_i 
            i_tranvec(3, i) = tranvec(3, i) - z_i 
          end do

          z_loop4: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop4: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop4: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = &
                            ior(img_j,ishft(cur_tran, 27))
#ifdef MPI
                          if (used_img_map(img_j) .eq. 0) then
                            used_img_map(img_j) = 1
                            used_img_cnt = used_img_cnt + 1
                            used_img_lst(used_img_cnt) = img_j
                          end if
#endif /* MPI */
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = &
                          ior(img_j,ishft(cur_tran,27))
#ifdef MPI
                        if (used_img_map(img_j) .eq. 0) then
                          used_img_map(img_j) = 1
                          used_img_cnt = used_img_cnt + 1
                          used_img_lst(used_img_cnt) = img_j
                        end if
#endif /* MPI */
                      end if
#ifdef MPI
                    else
                      if (used_img_ips_map(img_j) .eq. 0) then
                        used_img_ips_map(img_j) = 1
                        used_img_ips_cnt = used_img_ips_cnt + 1
                        used_img_ips_lst(used_img_ips_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do

                if (cur_bkt .eq. i_bkt) exit z_loop4 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop4
            end do y_loop4
          end do z_loop4
        end if
      end if
#else  /* NOT DIRFRC_COMTRANS */
      if (cutoffs_equal) then

        do i = 0, 17
          i_tranvec(1, i) = tranvec(1, i) - x_i 
          i_tranvec(2, i) = tranvec(2, i) - y_i 
          i_tranvec(3, i) = tranvec(3, i) - z_i 
        end do

        z_loop2: &
        do z_bkts_idx = z_bkts_lo, z_bkts_hi
          z_bkt = z_bkts(z_bkts_idx)
          z_tran = z_trans(z_bkts_idx)
          y_loop2: &
          do y_bkts_idx = y_bkts_lo, y_bkts_hi
            yz_bkt = z_bkt + y_bkts(y_bkts_idx)
            yz_tran = z_tran + y_trans(y_bkts_idx)
            x_loop2: &
            do x_bkts_idx = x_bkts_lo, x_bkts_hi
              cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
              cur_tran = x_trans(x_bkts_idx) + yz_tran
              do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                  if (excl_img_flags(img_j) .eq. 0) then
                    enc_img = ior(img_j, ishft(cur_tran, 27))
                    if (ico(iaci + img_iac(img_j)) .eq. 0) then
                      ee_eval_cnt = ee_eval_cnt + 1
                      img_j_ee_eval(ee_eval_cnt) = enc_img
                    else
                      full_eval_cnt = full_eval_cnt + 1
                      img_j_full_eval(full_eval_cnt) = enc_img
                    end if
#ifdef MPI
                    if (used_img_map(img_j) .eq. 0) then
                      used_img_map(img_j) = 1
                      used_img_cnt = used_img_cnt + 1
                      used_img_lst(used_img_cnt) = img_j
                    end if
                  else
                    if (used_img_ips_map(img_j) .eq. 0) then
                      used_img_ips_map(img_j) = 1
                      used_img_ips_cnt = used_img_ips_cnt + 1
                      used_img_ips_lst(used_img_ips_cnt) = img_j
                    end if
#endif /* MPI */
                  end if
                end if
              end do
              if (cur_bkt .eq. i_bkt) exit z_loop2 ! Done with cit buckets
                                                   ! image i claims pairs from.
            end do x_loop2
          end do y_loop2
        end do z_loop2

      else      ! cutoffs not equal

        do i = 0, 17
          i_tranvec(1, i) = tranvec(1, i) - x_i 
          i_tranvec(2, i) = tranvec(2, i) - y_i 
          i_tranvec(3, i) = tranvec(3, i) - z_i 
        end do

        z_loop4: &
        do z_bkts_idx = z_bkts_lo, z_bkts_hi
          z_bkt = z_bkts(z_bkts_idx)
          z_tran = z_trans(z_bkts_idx)
          y_loop4: &
          do y_bkts_idx = y_bkts_lo, y_bkts_hi
            yz_bkt = z_bkt + y_bkts(y_bkts_idx)
            yz_tran = z_tran + y_trans(y_bkts_idx)
            x_loop4: &
            do x_bkts_idx = x_bkts_lo, x_bkts_hi
              cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
              cur_tran = x_trans(x_bkts_idx) + yz_tran
              do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                r_sq = dx * dx + dy * dy + dz * dz
                if (r_sq .lt. cutlist_sq) then
                  if (excl_img_flags(img_j) .eq. 0) then
                    if (ico(iaci + img_iac(img_j)) .eq. 0) then
                      if (r_sq .lt. es_cut_sq) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = &
                          ior(img_j,ishft(cur_tran, 27))
#ifdef MPI
                        if (used_img_map(img_j) .eq. 0) then
                          used_img_map(img_j) = 1
                          used_img_cnt = used_img_cnt + 1
                          used_img_lst(used_img_cnt) = img_j
                        end if
#endif /* MPI */
                      end if
                    else
                      full_eval_cnt = full_eval_cnt + 1
                      img_j_full_eval(full_eval_cnt) = &
                        ior(img_j,ishft(cur_tran,27))
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
#ifdef MPI
                  else
                    if (used_img_ips_map(img_j) .eq. 0) then
                      used_img_ips_map(img_j) = 1
                      used_img_ips_cnt = used_img_ips_cnt + 1
                      used_img_ips_lst(used_img_ips_cnt) = img_j
                    end if
#endif /* MPI */
                  end if
                end if
              end do
              if (cur_bkt .eq. i_bkt) exit z_loop4 ! Done with cit buckets
                                                   ! image i claims pairs from.
            end do x_loop4
          end do y_loop4
        end do z_loop4

      end if
#endif /* DIRFRC_COMTRANS */

      total_pairs = total_pairs + full_eval_cnt + ee_eval_cnt       ! DBG

      if (dont_skip_belly_pairs) then
        call pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                          full_eval_cnt, img_j_full_eval, &
                          gbl_ipairs, num_packed)
      else
        call pack_nb_list_skip_belly_pairs(ee_eval_cnt, img_j_ee_eval, &
                                           full_eval_cnt, img_j_full_eval, &
                                           gbl_ipairs, igroup, img_atm_map, &
                                           num_packed)
      end if

      ! Clear excluded j images flags:

      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 0
      end do

      if (ifail .eq. 1) then
        flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi
        return
      end if

    end do ! end of images in i_bkt

    flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi

  end do ! end of all cit buckets with i images owned by this task

! write(0, *) 'DBG: total_pairs, listtot = ', &
!             total_pairs, num_packed
! write(0, *) 'DBG: avg packing = ', &
!             dble(total_pairs)/dble(num_packed)

  if (verbose .gt. 0) then
    if (master) write(mdout, *) 'listtot = ', num_packed
  end if

  return

contains

#include "pack_nb_list_dflt.i"

end subroutine get_nb_ips_list

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  check_my_atom_movement
!
! Description:  This routine checks if any atom on this task's atom list has
!               moved more than half the skin distance.  This is intended for
!               use only with parallel MD.
!
!*******************************************************************************

subroutine check_my_atom_movement(crd, saved_crd, my_atm_lst, skinnb, ntp, &
                                  new_list)

  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  double precision      :: crd(3, *)
  double precision      :: saved_crd(3, *)
  integer               :: my_atm_lst(*)
  double precision      :: skinnb
  integer               :: ntp
  logical               :: new_list

! Local variables:

  integer               :: atm_lst_idx
  double precision      :: dx, dy, dz
  double precision      :: box_dx, box_dy, box_dz
  double precision      :: distance_limit
  integer               :: n

! If anybody moves more than half the nbskin added to cutoff in generating the
! pairlist, update the pairlist.

  distance_limit = 0.25d0 * skinnb * skinnb

  if (ntp .eq. 0) then  ! Constant volume

    do atm_lst_idx = 1, my_atm_cnt
      n = my_atm_lst(atm_lst_idx)

      dx = crd(1, n) - saved_crd(1, n)
      dy = crd(2, n) - saved_crd(2, n)
      dz = crd(3, n) - saved_crd(3, n)

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  else  ! Constant pressure scaling.

    box_dx = pbc_box(1) / gbl_saved_box(1)
    box_dy = pbc_box(2) / gbl_saved_box(2)
    box_dz = pbc_box(3) / gbl_saved_box(3)

    do atm_lst_idx = 1, my_atm_cnt
      n = my_atm_lst(atm_lst_idx)

      dx = crd(1, n) - saved_crd(1, n) * box_dx
      dy = crd(2, n) - saved_crd(2, n) * box_dy
      dz = crd(3, n) - saved_crd(3, n) * box_dz

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  end if

  new_list = .false.

  return

end subroutine check_my_atom_movement
#endif /* MPI */

!*******************************************************************************
!
! Subroutine:  check_all_atom_movement
!
! Description:  This routine checks if any atom has moved more than half the
!               skin distance.  This is intended for all minimizations and
!               for single processor MD.
!*******************************************************************************

subroutine check_all_atom_movement(atm_cnt, crd, saved_crd, skinnb, ntp, &
                                   new_list)

  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, *)
  double precision      :: saved_crd(3, *)
  double precision      :: skinnb
  integer               :: ntp
  logical               :: new_list

! Local variables:

  double precision      :: dx, dy, dz
  double precision      :: box_dx, box_dy, box_dz
  double precision      :: distance_limit
  integer               :: n

! If anybody moves more than half the nbskin added to cutoff in generating the
! pairlist, update the pairlist.

  distance_limit = 0.25d0 * skinnb * skinnb

  if (ntp .eq. 0) then  ! Constant volume

    do n = 1, atm_cnt

      dx = crd(1, n) - saved_crd(1, n)
      dy = crd(2, n) - saved_crd(2, n)
      dz = crd(3, n) - saved_crd(3, n)

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  else  ! Constant pressure scaling.

    box_dx = pbc_box(1) / gbl_saved_box(1)
    box_dy = pbc_box(2) / gbl_saved_box(2)
    box_dz = pbc_box(3) / gbl_saved_box(3)

    do n = 1, atm_cnt

      dx = crd(1, n) - saved_crd(1, n) * box_dx
      dy = crd(2, n) - saved_crd(2, n) * box_dy
      dz = crd(3, n) - saved_crd(3, n) * box_dz

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  end if

  new_list = .false.
  return

end subroutine check_all_atom_movement

!*******************************************************************************
!
! Subroutine:  save_all_atom_crds
!
! Description:  This is needed for the skin test for buffered pairlists.
!*******************************************************************************

subroutine save_all_atom_crds(atm_cnt, crd, saved_crd)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: saved_crd(3, atm_cnt)

! Local variables:

  integer               :: i

  do i = 1, atm_cnt
    saved_crd(1, i) = crd(1, i)
    saved_crd(2, i) = crd(2, i)
    saved_crd(3, i) = crd(3, i)
  end do

  return

end subroutine save_all_atom_crds

!*******************************************************************************
!
! Subroutine:  save_imgcrds
!
! Description:  <TBS> 
!
!*******************************************************************************

#ifdef MPI
subroutine save_imgcrds(img_crd, saved_imgcrd, used_img_cnt, used_img_lst)
#else
subroutine save_imgcrds(img_cnt, img_crd, saved_imgcrd)
#endif

  use img_mod

  implicit none

! Formal arguments:

#ifdef MPI
  double precision      :: img_crd(3, *)
  double precision      :: saved_imgcrd(3, *)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
#else
  integer               :: img_cnt
  double precision      :: img_crd(3, *)
  double precision      :: saved_imgcrd(3, *)
#endif /* MPI */

! Local variables:

#ifdef MPI
  integer               :: lst_idx
#endif /* MPI */
  integer               :: img_idx

#ifdef MPI
  do lst_idx = 1, used_img_cnt
    img_idx = used_img_lst(lst_idx)
    saved_imgcrd(:, img_idx) = img_crd(:, img_idx)
  end do
#else
  do img_idx = 1, img_cnt
    saved_imgcrd(:,img_idx) = img_crd(:,img_idx)
  end do
#endif

  return

end subroutine save_imgcrds

!*******************************************************************************
!
! Subroutine:  adjust_imgcrds
!
! Description:  <TBS>
!              
!*******************************************************************************

#ifdef MPI
subroutine adjust_imgcrds(img_cnt, img_crd, img_atm_map, &
                          used_img_cnt, used_img_lst, &
                          saved_imgcrd, crd, saved_crd, ntp)
#else
subroutine adjust_imgcrds(img_cnt, img_crd, img_atm_map, &
                          saved_imgcrd, crd, saved_crd, ntp)
#endif /* MPI */

  use gbl_datatypes_mod
  use img_mod
  use pbc_mod

  implicit none

! Formal arguments:

#ifdef MPI
  integer               :: img_cnt
  double precision      :: img_crd(3, *)
  integer               :: img_atm_map(*)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
  double precision      :: saved_imgcrd(3, *)
  double precision      :: crd(3, *)
  double precision      :: saved_crd(3, *)
  integer               :: ntp
#else
  integer               :: img_cnt
  double precision      :: img_crd(3, *)
  integer               :: img_atm_map(*)
  double precision      :: saved_imgcrd(3, *)
  double precision      :: crd(3, *)
  double precision      :: saved_crd(3, *)
  integer               :: ntp
#endif /* MPI */

! Local variables:

  integer               :: atm_idx
  integer               :: img_idx
  double precision      :: box_del(3)
#ifdef MPI
  integer               :: lst_idx
#endif /* MPI */

#ifdef MPI
  if (ntp .eq. 0) then  ! Constant volume

      do lst_idx = 1, used_img_cnt
        img_idx = used_img_lst(lst_idx)
        atm_idx = img_atm_map(img_idx)
        img_crd(:,img_idx) = crd(:,atm_idx) + saved_imgcrd(:,img_idx) - &
                             saved_crd(:,atm_idx)
      end do

  else  ! Constant pressure scaling.

    box_del(:) = pbc_box(:) / gbl_saved_box(:)

      do lst_idx = 1, used_img_cnt
        img_idx = used_img_lst(lst_idx)
        atm_idx = img_atm_map(img_idx)
        img_crd(:,img_idx) = crd(:,atm_idx) + (saved_imgcrd(:,img_idx) - &
                             saved_crd(:,atm_idx)) * box_del(:)
      end do

  end if

#else

  if (ntp .eq. 0) then  ! Constant volume

    do img_idx = 1, img_cnt
        atm_idx = img_atm_map(img_idx)
        img_crd(:,img_idx) = crd(:,atm_idx) + saved_imgcrd(:,img_idx) - &
                             saved_crd(:,atm_idx)
    end do

  else  ! Constant pressure scaling.

    box_del(:) = pbc_box(:) / gbl_saved_box(:)

    do img_idx = 1, img_cnt
      atm_idx = img_atm_map(img_idx)
      img_crd(:,img_idx) = crd(:,atm_idx) + (saved_imgcrd(:,img_idx) - &
                           saved_crd(:,atm_idx)) * box_del(:)
    end do

  end if

#endif /* MPI */

  return

end subroutine adjust_imgcrds

end module nb_pairlist_mod
