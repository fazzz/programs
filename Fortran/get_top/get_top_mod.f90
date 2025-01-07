! module to read the topology file of gromacs

module get_gmx_top
  implicit none

  type GMX_data
! default
!     integer nbfunc, comb_rule
!     character(4) gen_pairs
!     double precision fudgeLJ, fudgeQQ

! atomtypes
!     character,allocatable,dimension(:,3) :: atomtype_name, bondtype
!     double precision,allocatable,dimension(:) :: atomtype_mass, atomtype_charge
!     character,allocatable,dimension(:) :: ptye
!     double precision,allocatable,dimension(:) :: atomtype_sigma, atomtype_epsilon

! atoms
!     integer,allocatable,dimension(:) :: nrexcl
     character,allocatable,dimension(:,3) :: atom_type
     integer,allocatable,dimension(:) :: atom_resi
     character,allocatable,dimension(;,3) :: atom_res
     character,allocatable,dimension(;,4) :: atom_name
     double precision,allocatable,dimension(:) :: atom_charge
     double precision,allocatable,dimension(:) :: atom_mass

! position_restraint
!     integer,allocatable,dimension(:) :: pr_i, pr_funct
!     double precision,allocatable,dimension(:) :: fcx, fcy, fcz

! system
     character,allocatable,dimension(:,100) :: system_tytle

! molecules
     character,allocatable,dimension(:,10) :: Compound
     integer nmols
     
! molecule_type
!     character,allocatable,dimension(:,10) :: name
!     integer,allocatable,dimension(:) :: nrexcl

! bonds
     integer,allocatable,dimension(:) :: bi, bj, b_funct
     double precision,allocatable,dimension(:) :: br
     integer,allocatable,dimension(:) :: bk

! pairs
     integer,allocatable,dimension(:) :: pi, pj
     integer,allocatable,dimension(:) :: p_funct

! angles
     integer,allocatable,dimension(:) :: ai, aj, ak
     integer,allocatable,dimension(:) :: a_funct
     double precision,allocatable,dimension(:) :: theta, cth

! dihedrals
     integer,allocatable,dimension(:) :: di, dj, dk, dl
     integer,allocatable,dimension(:) :: d_funct
     double precision,allocatable,dimension(:) :: c0, c1, c2, c3, c4

  end type GMX_data

contains

  subroutine get_top( file_top, num_of_residue, num_of_atom, GMX_data )
    implicit none
!
    character(len=200),intent(in) :: file_top
    integer :: num_of_atom, num_of_residue
    type (GMX_data),intent(inout) :: GMX_data

    real(8) :: rtemp1, rtemp2
    real(8) :: charge(num_of_atom), mass(num_of_atom)
    character(len=3) :: ctemp3, resname(num_of_residue)
    character(len=4) :: ctemp4
    character(len=4) :: atomname(num_of_atom)
    character :: FLAG
!
    open(10,file=file_top)
    do
       if (MODE .eq. A) then
          read(10,("A1")),FLAG
          if (FLAG .eq. "#") then
             
          endif
          if (FLAG .eq. "[") then
          
          endif
       else if (MODE .eq. B) then

       else then
       endif

    enddo
  end subroutine get_top
end module get_gmx_top
