! This module provides the subroutine reallocate 
! which is used to re-size the array.
module reallocate_Fortran
  implicit none

  contains

    subroutine reallocate_real(A, newsize)
      implicit none
      real, allocatable, intent(inout) :: A
      integer, intent(in) :: newsize

      real, allocatable, intent(inout) :: B

      allocate(B(LBOUND(A):UBOUND(A)))
      
      B = A

      deallocate(A)

      allocate(A(newsize))

      A(LBOUND(B):UBOUND(B)) = B

      deallocate(B)
    end subroutine reallocate_real
end module reallocate_Fortran
