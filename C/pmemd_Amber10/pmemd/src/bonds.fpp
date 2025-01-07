#include "copyright.i"

!*******************************************************************************
!
! Module:  bonds_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module bonds_mod

  use gbl_datatypes_mod

  implicit none

! The following are derived from prmtop bond info:

  integer, save                         :: cit_nbonh, cit_nbona

  type(bond_rec), allocatable, save     :: cit_h_bond(:)
  type(bond_rec), allocatable, save     :: cit_a_bond(:)

contains

!*******************************************************************************
!
! Subroutine:  bonds_setup
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine bonds_setup(num_ints, num_reals, use_atm_map)

  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals
  integer                       :: use_atm_map(natom)

! Local variables:

  integer               :: alloc_failed
  type(bond_rec)        :: bonds_copy(nbonh + nbona)

  ! This routine can handle reallocation, and thus can be called multiple
  ! times.

  call find_my_bonds(nbonh, gbl_bond, cit_nbonh, bonds_copy, use_atm_map)

  if (cit_nbonh .gt. 0) then
    if (allocated(cit_h_bond)) then
      if (size(cit_h_bond) .lt. cit_nbonh) then
        num_ints = num_ints - size(cit_h_bond) * bond_rec_ints
        deallocate(cit_h_bond)
        allocate(cit_h_bond(cit_nbonh), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_h_bond) * bond_rec_ints
      end if
    else
      allocate(cit_h_bond(cit_nbonh), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_h_bond) * bond_rec_ints
    end if
    cit_h_bond(1:cit_nbonh) = bonds_copy(1:cit_nbonh)
  end if

  call find_my_bonds(nbona, gbl_bond(bonda_idx), cit_nbona, bonds_copy, &
                         use_atm_map)

  if (cit_nbona .gt. 0) then
    if (allocated(cit_a_bond)) then
      if (size(cit_a_bond) .lt. cit_nbona) then
        num_ints = num_ints - size(cit_a_bond) * bond_rec_ints
        deallocate(cit_a_bond)
        allocate(cit_a_bond(cit_nbona), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_a_bond) * bond_rec_ints
      end if
    else
      allocate(cit_a_bond(cit_nbona), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_a_bond) * bond_rec_ints
    end if
    cit_a_bond(1:cit_nbona) = bonds_copy(1:cit_nbona)
  end if

  return

end subroutine bonds_setup

!*******************************************************************************
!
! Subroutine:  find_my_bonds
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine find_my_bonds(bond_cnt, bonds, my_bond_cnt, my_bonds, &
                             use_atm_map)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: bond_cnt
  type(bond_rec)        :: bonds(bond_cnt)
  integer               :: my_bond_cnt
  type(bond_rec)        :: my_bonds(*)
  integer               :: use_atm_map(*)

! Local variables:

  integer               :: atm_i, atm_j, bonds_idx

! Find all bonds for which this process owns either atom:

  my_bond_cnt = 0

  do bonds_idx = 1, bond_cnt

    atm_i = bonds(bonds_idx)%atm_i
    atm_j = bonds(bonds_idx)%atm_j

#ifdef MPI
    if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
      my_bond_cnt = my_bond_cnt + 1
      my_bonds(my_bond_cnt) = bonds(bonds_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
    else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
      my_bond_cnt = my_bond_cnt + 1
      my_bonds(my_bond_cnt) = bonds(bonds_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
    end if
#else
    my_bond_cnt = my_bond_cnt + 1
    my_bonds(my_bond_cnt) = bonds(bonds_idx)
    use_atm_map(atm_i) = 1
    use_atm_map(atm_j) = 1
#endif

  end do

  return

end subroutine find_my_bonds

!*******************************************************************************
!
! Subroutine:  get_bond_energy
!
! Description:
!              
! Routine to get bond energy and forces for the potential of cb*(b-b0)**2.
!
!*******************************************************************************

subroutine get_bond_energy(bond_cnt, bond, x, frc, bond_energy)

  use nmr_calls_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: bond_cnt
  type(bond_rec)        :: bond(*)
  double precision      :: x(3, *)
  double precision      :: frc(3, *)
  double precision      :: bond_energy

! Local variables:

  double precision      :: da
  double precision      :: df
  double precision      :: dfw
  integer               :: i, j, ic, jn
  double precision      :: lcl_bond_energy
  double precision      :: xa, ya, za
  double precision      :: rij
  double precision      :: xij, yij, zij

  lcl_bond_energy = 0.0d0

! Grand loop for the bond stuff:

  do jn = 1, bond_cnt

    i = bond(jn)%atm_i
    j = bond(jn)%atm_j

! Calculation of the bond vector:

    xij = x(1, i) - x(1, j)
    yij = x(2, i) - x(2, j)
    zij = x(3, i) - x(3, j)

    rij = sqrt(xij * xij + yij * yij + zij * zij)

! Calculation of the energy and deriv:

    ic = bond(jn)%parm_idx
    da = rij - gbl_req(ic)

#ifdef MPI
! If ebdev is ever supported under mpi, you will need to treat it like the
! energy...
#else
    ebdev = ebdev + da * da ! For rms deviation from ideal bonds:
#endif
    
    df = gbl_rk(ic) * da
    dfw = (df + df) / rij

! Calculation of the force:

    xa = dfw * xij
    ya = dfw * yij
    za = dfw * zij

#ifdef MPI
    ! We use atm_i to determine who sums up the energy...
    if (gbl_atm_owner_map(i) .eq. mytaskid) then
      lcl_bond_energy = lcl_bond_energy + df * da
      frc(1, i) = frc(1, i) - xa
      frc(2, i) = frc(2, i) - ya
      frc(3, i) = frc(3, i) - za
    end if

    if (gbl_atm_owner_map(j) .eq. mytaskid) then
      frc(1, j) = frc(1, j) + xa
      frc(2, j) = frc(2, j) + ya
      frc(3, j) = frc(3, j) + za
    end if
#else
      lcl_bond_energy = lcl_bond_energy + df * da
      frc(1, i) = frc(1, i) - xa
      frc(2, i) = frc(2, i) - ya
      frc(3, i) = frc(3, i) - za
      frc(1, j) = frc(1, j) + xa
      frc(2, j) = frc(2, j) + ya
      frc(3, j) = frc(3, j) + za
#endif

  end do

  bond_energy = lcl_bond_energy

  return

end subroutine get_bond_energy

end module bonds_mod
