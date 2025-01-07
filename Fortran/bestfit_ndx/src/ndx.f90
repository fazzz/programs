module ndx
  implicit none

contains

  subroutine read_ndx (unitnum, Nndx, atom_list)

    implicit none
    integer, intent(in) :: unitnum
    integer, intent(inout) :: Nndx
    integer, allocatable, dimension(:), intent(out) :: atom_list

    integer :: i, j, k
    integer, allocatable, dimension(:) :: B
    integer v, err
    
    i = 0
    read (unitnum, '()')
    do 
       read (unitnum, '(I4,1X)', advance="NO",iostat=err, end=999), v

       if (err .eq. 0 ) then
          if ( i > 0) then
             allocate(B(i))
             B = atom_list
             deallocate(atom_list)
          end if
          i = i + 1
          allocate(atom_list(i))
          if ( i - 1 > 0) then
             atom_list(1:i-1) = B(1:i-1)
             deallocate(B)
          end if
          atom_list(i) = v
       end if
    end do
999 Nndx = i
    
!  do i = 1,Nndx,1
!     write(*,'(I4,1X)', advance="NO"), atom_list(i)
!     if (mod(i,10) .eq. 0) then
!        write(*,"()")
!     end if
!  end do

  end subroutine read_ndx
      
end module ndx
