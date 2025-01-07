! This modeule contains the function to read
! topology files (gromacs)

module gettop
implicit none

contains
  subroutine get_top( unitnum, NA, nr,  atype,  resid,  res,  atom,  cgnr, charg,  mass )
    implicit none

    integer,intent(in) :: unitnum
    integer,intent(in) :: NA
    integer,allocatable,dimension(:),intent(inout) :: nr, resid, cgnr
    double precision,allocatable,dimension(:),intent(inout) :: charg, mass
!    character,allocatable,dimension(:,:),intent(inout) :: atom
    character(5),intent(inout) :: atype, res, atom

    integer i

    write(*,*) NA

    do i=1,NA
       read(unitnum,'(I6,A5,I6,1x,A5,1x,A5,I5,4x,F9.6,5x,F8.5)'),&
            nr(i), atype, resid(i), res, atom, cgnr(i), charg(i), mass(i)
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! write(*,*), &                                                      ! !
       !      nr(i), atype, resid(i), res, atom, cgnr(i), charg(i), mass(i) ! !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
    enddo

  end subroutine get_top
end module gettop
