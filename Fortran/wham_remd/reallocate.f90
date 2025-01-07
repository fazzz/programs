! This module provides the subroutine reallocate 
! which is used to re-size the array.
module reallocate
  implicit none

  contains

    subroutine reallocate_real(A, newsize)
      implicit none
      real, allocatable, dimension(:), intent(inout) :: A
      integer,  intent(in) :: newsize

      real, allocatable, dimension(:) :: B

      allocate(B(lbound(A,1):ubound(A,1)))
      
      B = A

      deallocate(A)

      allocate(A(newsize))

      A(lbound(B,1):ubound(B,1)) = B

      deallocate(B)
    end subroutine reallocate_real
end module reallocate
