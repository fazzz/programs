module ndx
  implicit none

contains

  subroutine read_ndx (unitnum, Nndx, atom_list)

    implicit none
    integer, intent(in) :: unitnum
    integer, intent(inout) :: Nndx
    integer, allocatable, dimension(:), intent(inout) :: atom_list

    integer :: i, j, k
    integer, allocatable, dimension(:) :: B
    integer v(15)

    write(*,'("YES")')
    
    i = 1
    read (unitnum, '()')
    do 
       read (unitnum,'(I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4,1X,I4)',&
            end=999), v(1),v(2),v(3),v(4),v(5),v(11),v(12),v(13),v(14),v(15)
       do k = 1,15,1
          allocate(B(i))
          B = atom_list
          i = i + 1
          deallocate(atom_list)
          allocate(atom_list(i))
          atom_list(1:i-1) = B(i-1)
          deallocate(B)
          atom_list(i) = v(k)
       end do
    end do
999 Nndx = i - 1
    
  end subroutine read_ndx
      
end module ndx
