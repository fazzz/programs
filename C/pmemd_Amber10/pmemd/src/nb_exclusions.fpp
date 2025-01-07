#include "copyright.i"

!*******************************************************************************
!
! Module: nb_exclusions_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module nb_exclusions_mod

  use gbl_datatypes_mod

  implicit none

  type(listdata_rec), allocatable, save :: atm_nb_maskdata(:)
  integer, allocatable, save            :: atm_nb_mask(:)
  integer, allocatable, save            :: gbl_nb_adjust_pairlst(:)

contains

!*******************************************************************************
!
! Subroutine:  alloc_nb_exclusions_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_nb_exclusions_mem(atm_cnt, ext_cnt, num_ints, num_reals)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt
  integer                       :: ext_cnt

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed

  allocate(gbl_nb_adjust_pairlst(atm_cnt + ext_cnt + ext_cnt), &
           atm_nb_maskdata(atm_cnt), &
           atm_nb_mask(ext_cnt + ext_cnt), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(gbl_nb_adjust_pairlst) + &
                        2 * size(atm_nb_maskdata) + &
                        size(atm_nb_mask)

  gbl_nb_adjust_pairlst(:) = 0
  atm_nb_maskdata(:) = listdata_rec(0,0)
  atm_nb_mask(:) = 0

  return

end subroutine alloc_nb_exclusions_mem

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_nb_exclusions_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_nb_exclusions_dat(atm_cnt, ext_cnt)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer       :: atm_cnt
  integer       :: ext_cnt

! Local variables:

  integer       :: num_ints, num_reals  ! returned values discarded

  ! Nothing to broadcast.  We just allocate storage in the non-master nodes.

  if (.not. master) then
    num_ints = 0
    num_reals = 0
    call alloc_nb_exclusions_mem(atm_cnt, ext_cnt, num_ints, num_reals)
  end if

  ! The allocated data is not initialized from the master node.

  return

end subroutine bcast_nb_exclusions_dat
#endif

!*******************************************************************************
!
! Subroutine:  make_atm_excl_mask_list
!
! Description:  This routine takes as input the numex() and natex() arrays from
!               the prmtop.  The prmtop format basically consists of a count of
!               excluded nonbond interactions for each atom (in numex()) and a
!               correlated concatenated "list of sublists" of the actual
!               exclusions in natex() (see the prmtop format doc to understand
!               the input - but the doc is not completely correct :-( ))).
!               The input format is flawed in that it is one-directional (only
!               listing exclusions for atoms with a higher atom number) and also
!               it is not possible to get data for an individual atom without
!               traversing the numex() array.  Here we basically double the list
!               (listing nonbonded interactions in both directions) and also
!               change the list data to allow finding the data for a given atom
!               in one probe (the listdata_rec struct is used, and we make an
!               atm_nb_maskdata() array of this structure and an atm_nb_mask()
!               array that is equivalent to natex() but that accounts for
!               nonbonded interactions in both directions).
!              
!*******************************************************************************

subroutine make_atm_excl_mask_list(atm_cnt, numex, next, natex)

  use gbl_constants_mod
  use file_io_dat_mod, only : mdout
  use constraints_mod, only : atm_igroup
  use mdin_ctrl_dat_mod, only : ibelly
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  integer               :: numex(atm_cnt) ! Excluded atom count for each atom.
  integer               :: next           ! natex() size, from prmtop.
  integer               :: natex(next)    ! Excluded atom concatenated list.

! Local variables:

  integer               :: atm_i, atm_j
  integer               :: lst_idx, sublst_idx, num_sublst
  integer               :: mask_idx
  integer               :: offset
  integer               :: total_excl

! Double the mask to deal with our list generator

! Pass 1: get pointers, check size

  lst_idx = 0

  atm_nb_maskdata(:)%cnt = 0       ! array assignment

  do atm_i = 1, atm_cnt - 1         ! last atom never has any...
    num_sublst = numex(atm_i)
    do sublst_idx = 1, num_sublst
      atm_j = natex(lst_idx + sublst_idx)
      if (atm_j .gt. 0) then
        atm_nb_maskdata(atm_i)%cnt = atm_nb_maskdata(atm_i)%cnt + 1
        atm_nb_maskdata(atm_j)%cnt = atm_nb_maskdata(atm_j)%cnt + 1
      end if
    end do
    lst_idx = lst_idx + num_sublst
  end do

  total_excl = 0

  do atm_i = 1, atm_cnt
    total_excl = total_excl + atm_nb_maskdata(atm_i)%cnt
  end do

  if (total_excl .gt. 2 * next) then
    write(mdout, '(a,a)') error_hdr, &
      'Variable next has inconsistent value in prmtop!'
    call mexit(6, 1)
  end if

  offset = 0

  do atm_i = 1, atm_cnt
    atm_nb_maskdata(atm_i)%offset = offset
    offset = offset + atm_nb_maskdata(atm_i)%cnt
  end do

! Pass 2: fill mask array

  lst_idx = 0

  atm_nb_maskdata(:)%cnt = 0       ! array assignment
  
  if (ibelly .eq. 0) then
    do atm_i = 1, atm_cnt - 1
      num_sublst = numex(atm_i)
      do sublst_idx = 1, num_sublst
        atm_j = natex(lst_idx + sublst_idx)
        if (atm_j .gt. 0) then

          atm_nb_maskdata(atm_j)%cnt = atm_nb_maskdata(atm_j)%cnt + 1
          mask_idx = atm_nb_maskdata(atm_j)%offset + &
                     atm_nb_maskdata(atm_j)%cnt
          atm_nb_mask(mask_idx) = atm_i

          atm_nb_maskdata(atm_i)%cnt = atm_nb_maskdata(atm_i)%cnt + 1
          mask_idx = atm_nb_maskdata(atm_i)%offset + &
                     atm_nb_maskdata(atm_i)%cnt
          atm_nb_mask(mask_idx) = atm_j

        end if
      end do
      lst_idx = lst_idx + num_sublst
    end do
  else
    do atm_i = 1, atm_cnt - 1
      num_sublst = numex(atm_i)
      do sublst_idx = 1, num_sublst
        atm_j = natex(lst_idx + sublst_idx)
        if (atm_j .gt. 0) then
          if (atm_igroup(atm_i) .ne. 0 .and. atm_igroup(atm_j) .ne. 0) then

            atm_nb_maskdata(atm_j)%cnt = atm_nb_maskdata(atm_j)%cnt + 1
            mask_idx = atm_nb_maskdata(atm_j)%offset + &
                       atm_nb_maskdata(atm_j)%cnt
            atm_nb_mask(mask_idx) = atm_i

            atm_nb_maskdata(atm_i)%cnt = atm_nb_maskdata(atm_i)%cnt + 1
            mask_idx = atm_nb_maskdata(atm_i)%offset + &
            atm_nb_maskdata(atm_i)%cnt
            atm_nb_mask(mask_idx) = atm_j

          end if
        end if
      end do
      lst_idx = lst_idx + num_sublst
    end do
  end if

  return

end subroutine make_atm_excl_mask_list

!*******************************************************************************
!
! Subroutine:   make_nb_adjust_pairlst
!
! Description:  Make a list of pairs for which we must make a reciprocal space
!               nonbonded adjustment.  We here are creating a list that is
!               used similarly to the bond, angle, dihedral lists.
!
!*******************************************************************************

#ifdef MPI
subroutine make_nb_adjust_pairlst(my_atm_cnt, my_atm_lst, &
                                  atm_owner_map, use_atm_map, &
                                  nb_adjust_pairlst, atm_maskdata, atm_mask)
#else
subroutine make_nb_adjust_pairlst(atm_cnt, use_atm_map, &
                                  nb_adjust_pairlst, atm_maskdata, atm_mask)
#endif

  use gbl_datatypes_mod
  use img_mod
#ifdef MPI
  use parallel_dat_mod, only : mytaskid
#endif
  use pbc_mod

  implicit none

! Formal arguments:

#ifdef MPI
  integer               :: my_atm_cnt
  integer               :: my_atm_lst(my_atm_cnt)
  integer               :: atm_owner_map(*)
#else
  integer               :: atm_cnt
#endif
  integer               :: use_atm_map(*)
  integer               :: nb_adjust_pairlst(*)
  type(listdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)

! Local variables:

#ifdef MPI
  integer               :: atm_lst_idx
#endif
  integer               :: atm_i, atm_j
  integer               :: atm_mask_off, atm_mask_idx
  integer               :: excl_atm_j_cnt
  integer               :: nb_adj_pairlst_len

  nb_adj_pairlst_len = 0

  ! Loop over the atoms you own, finding excluded atoms you should process in
  ! nb_adjust().  The uniprocessor code sets use_atm_map(), presumably
  ! unnecessarily except perhaps for debugging; this was done here to be
  ! consistent with other code that sets use_atm_map().  This code should only
  ! execute at program setup and when the atoms are redistributed, which is
  ! an infrequent event.

#ifdef MPI
  ! Only atm_i owned by this task are processed.  For atm_j, if it is not owned
  ! then the atm id is negated as a flag to the nb_adjust() code that special
  ! processing is required.  In this case, the task will calc forces only for
  ! the atom it owns, and will only calc energies if atm_i < atm_j.  All atm_j
  ! not owned must be put on the extra_used_atoms list (or otherwise be
  ! claimed to insure that we keep the coordinates updated).
#endif /* MPI */

#ifdef MPI
  do atm_lst_idx = 1, my_atm_cnt
    atm_i = my_atm_lst(atm_lst_idx)
#else
  do atm_i = 1, atm_cnt
#endif
      use_atm_map(atm_i) = 1
      nb_adj_pairlst_len = nb_adj_pairlst_len + 1
      excl_atm_j_cnt = 0

      atm_mask_off = atm_maskdata(atm_i)%offset

      do atm_mask_idx = atm_mask_off + 1, &
                        atm_mask_off + atm_maskdata(atm_i)%cnt

        atm_j = atm_mask(atm_mask_idx)

#ifdef MPI
        if (atm_owner_map(atm_j) .eq. mytaskid) then
#endif
          if (atm_i .lt. atm_j) then
            excl_atm_j_cnt = excl_atm_j_cnt + 1
            nb_adjust_pairlst(nb_adj_pairlst_len + excl_atm_j_cnt) = atm_j
            use_atm_map(atm_j) = 1
          end if
#ifdef MPI
        else
          excl_atm_j_cnt = excl_atm_j_cnt + 1
          ! Note negation to mark this as a mixed ownership pair...
          nb_adjust_pairlst(nb_adj_pairlst_len + excl_atm_j_cnt) = -atm_j
          use_atm_map(atm_j) = 1
        end if
#endif

      end do

      nb_adjust_pairlst(nb_adj_pairlst_len) = excl_atm_j_cnt
      nb_adj_pairlst_len = nb_adj_pairlst_len + excl_atm_j_cnt

  end do ! main atom loop

  return

end subroutine make_nb_adjust_pairlst

!*******************************************************************************
!
! Subroutine:  nb_adjust
!
! Description:  The part of ewald due to gaussian counterion about an atom you
!               are bonded to or otherwise for which you do not compute the
!               nonbond pair force. NECESSARY since you ARE computing this pair
!               in the reciprocal sum.
!*******************************************************************************

#ifdef MPI
subroutine nb_adjust(my_atm_cnt, my_atm_lst, &
                     charge, crd, nb_adj_pairlst, eed_cub, &
                     frc, ene, virial)
#else
subroutine nb_adjust(atm_cnt, &
                     charge, crd, nb_adj_pairlst, eed_cub, &
                     frc, ene, virial)
#endif

  use mdin_ewald_dat_mod
#ifdef MPI
  use parallel_dat_mod, only : mytaskid
#endif

  implicit none

! Formal arguments:

#ifdef MPI
  integer               :: my_atm_cnt
  integer               :: my_atm_lst(my_atm_cnt)
#else
  integer               :: atm_cnt
#endif
  double precision      :: charge(*)
  double precision      :: crd(3, *)
  integer               :: nb_adj_pairlst(*)
  double precision      :: eed_cub(*)
  double precision      :: frc(3, *)
  double precision      :: ene
  double precision      :: virial(3, 3)

! Local variables:

  integer               :: atm_i, atm_j
#ifdef MPI
  integer               :: atm_lst_idx
  logical               :: owned_atm_j
#endif
  integer               :: excl_atm_j_cnt
  integer               :: nxt_sublst
  integer               :: sublst_idx
  double precision      :: ene_stk
  double precision      :: ew_coeff_stk
  double precision      :: eedtbdns_stk
  double precision      :: delx, dely, delz, delr, delr2, delrinv
  double precision      :: cgi, cgi_cgj
  double precision      :: del, dx, x
  double precision      :: erfc, derfc, d0, d1
  double precision      :: df, dfx, dfy, dfz
  double precision      :: vxx, vxy, vxz, vyy, vyz, vzz
  double precision      :: x_i, y_i, z_i
  integer               :: ind

  double precision, parameter   :: half = 1.d0/2.d0
  double precision, parameter   :: third = 1.d0/3.d0

#ifdef MPI
  if (my_atm_cnt .eq. 0) return
#endif

  ene_stk = 0.d0
  ew_coeff_stk = ew_coeff
  eedtbdns_stk = eedtbdns

  vxx = 0.d0
  vxy = 0.d0
  vxz = 0.d0
  vyy = 0.d0
  vyz = 0.d0
  vzz = 0.d0

  del = 1.d0 / eedtbdns

  nxt_sublst = 0

#ifdef MPI
  owned_atm_j = .true.

  do atm_lst_idx = 1, my_atm_cnt
    atm_i = my_atm_lst(atm_lst_idx)
#else
  do atm_i = 1, atm_cnt
#endif

    cgi = charge(atm_i)

    nxt_sublst = nxt_sublst + 1
    excl_atm_j_cnt = nb_adj_pairlst(nxt_sublst)

    x_i = crd(1, atm_i)
    y_i = crd(2, atm_i)
    z_i = crd(3, atm_i)

    do sublst_idx = nxt_sublst + 1, nxt_sublst + excl_atm_j_cnt

      atm_j = nb_adj_pairlst(sublst_idx)

#ifdef MPI
      if (atm_j .lt. 0) then
        owned_atm_j = .false.
        atm_j = -atm_j
      end if
#endif

      delx = crd(1, atm_j) - x_i
      dely = crd(2, atm_j) - y_i
      delz = crd(3, atm_j) - z_i

      delr2 = delx * delx + dely * dely + delz * delz

      ! Similar code to that in pairs calc code; however the only valid option
      ! here is the erfc switch, and dxdr = ewaldcof.

      delr = sqrt(delr2)
      delrinv = 1.d0 / delr

      x = ew_coeff_stk * delr
      ind = int(eedtbdns_stk * x)
      dx = x - dble(ind) * del
      ind = ishft(ind, 2)             ! 4 * ind
      
      ! Cubic spline on erfc derfc:

      erfc = eed_cub(1 + ind) + dx * (eed_cub(2 + ind) + &
             dx * (eed_cub(3 + ind) + dx * eed_cub(4 + ind) * third) * half)

      derfc = eed_cub(2 + ind) + dx * (eed_cub(3 + ind) + &
              dx * eed_cub(4 + ind) * half)

      d0 = (erfc - 1.d0) * delrinv
      d1 = (d0 - ew_coeff_stk * derfc) * delrinv * delrinv

      cgi_cgj = cgi * charge(atm_j)

      df = cgi_cgj * d1

      dfx = delx * df
      dfy = dely * df
      dfz = delz * df

      frc(1, atm_i) = frc(1, atm_i) - dfx
      frc(2, atm_i) = frc(2, atm_i) - dfy
      frc(3, atm_i) = frc(3, atm_i) - dfz

#ifdef MPI
      if (owned_atm_j) then
#endif

        frc(1, atm_j) = frc(1, atm_j) + dfx
        frc(2, atm_j) = frc(2, atm_j) + dfy
        frc(3, atm_j) = frc(3, atm_j) + dfz

        ene_stk = ene_stk + cgi_cgj * d0
        vxx = vxx - dfx * delx
        vxy = vxy - dfx * dely
        vxz = vxz - dfx * delz
        vyy = vyy - dfy * dely
        vyz = vyz - dfy * delz
        vzz = vzz - dfz * delz

#ifdef MPI
      else

        if (atm_i .lt. atm_j) then

          ene_stk = ene_stk + cgi_cgj * d0
          vxx = vxx - dfx * delx
          vxy = vxy - dfx * dely
          vxz = vxz - dfx * delz
          vyy = vyy - dfy * dely
          vyz = vyz - dfy * delz
          vzz = vzz - dfz * delz

        end if

        owned_atm_j = .true.    ! default setting...

      end if
#endif

    end do

    nxt_sublst = nxt_sublst + excl_atm_j_cnt

  end do

  virial(1, 1) = vxx
  virial(1, 2) = vxy
  virial(2, 1) = vxy
  virial(1, 3) = vxz
  virial(3, 1) = vxz
  virial(2, 2) = vyy
  virial(2, 3) = vyz
  virial(3, 2) = vyz
  virial(3, 3) = vzz

  ene = ene_stk

  return

end subroutine nb_adjust

end module nb_exclusions_mod
