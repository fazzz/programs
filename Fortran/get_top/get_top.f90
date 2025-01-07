! This modeule contains the function to read
! topology files (gromacs)
subroutine get_top( file_top, num_of_residue, &
     num_of_atom, atomname, charge, mass ,resname )
implicit none
!
integer :: i, iatm, num_of_atom, num_of_residue, ires
real(8) :: rtemp1, rtemp2
real(8) :: charge(num_of_atom), mass(num_of_atom)
character(len=3) :: ctemp3, resname(num_of_residue)
character(len=4) :: ctemp4
character(len=200) :: file_top
character(len=4) :: atomname(num_of_atom)
!
open(10,file=file_top)
do i = 1, 26
  read(10,*)
enddo
!
read(10,'(T15,A3)') ctemp3
ires = 1
resname(ires) = ctemp3
!
iatm = 0
do while ( iatm.ne.num_of_atom )
  read(10,100) ctemp3, ctemp4, rtemp1, rtemp2
  if ( ctemp4.ne."    " ) then
       iatm = iatm + 1
       atomname(iatm) = ctemp4
       charge(iatm) = rtemp1
       mass(iatm) = rtemp2
       if ( atomname(iatm).eq."   O" .or. atomname(iatm)(2:2).eq."O" ) then
            mass(iatm) = mass(iatm) * 100.0d0
       endif
  else
           ires = ires + 1
           resname(ires) = ctemp3
           print*, ires, resname(ires)
  endif
enddo
close(10)
100 format((T15,A3),(T35,A4),(T50,F7.4),(T63,F5.2))
!
stop
end subroutine
!
