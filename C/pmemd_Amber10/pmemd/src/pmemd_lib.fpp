#include "copyright.i"

!*******************************************************************************
!
! Module: pmemd_lib_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module pmemd_lib_mod

use file_io_dat_mod

  implicit none

#ifdef MPI
  ! This needs to be visible for adding groups; sigh...
  integer, save         :: world_group
#endif /* MPI */

contains

!*******************************************************************************
!
! Subroutine:  mstartup
!
! Description: <TBS>
!              
!*******************************************************************************

#ifdef MPI
subroutine mstartup(mytaskid, numtasks)

#ifdef USE_MPI_MODULE
  USE_MPI_MODULE
  implicit none
#else
  implicit none
#include <mpif.h>
#endif

! Formal arguments:

  integer, intent(out)  :: mytaskid
  integer, intent(out)  :: numtasks

! Local variables:

  integer       :: err_ret_code

  call mpi_init(err_ret_code)
  call mpi_comm_rank(mpi_comm_world, mytaskid, err_ret_code)
  call mpi_comm_size(mpi_comm_world, numtasks, err_ret_code)
  call mpi_comm_group(mpi_comm_world, world_group, err_ret_code)

  return

end subroutine mstartup
#endif /* MPI */

!*******************************************************************************
!
! Subroutine:  mexit
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine mexit(ifil, i)

#ifdef MPI
#ifdef USE_MPI_MODULE
  USE_MPI_MODULE
  implicit none
#else
  implicit none
#include <mpif.h>
#endif
#endif

! Formal arguments:

  integer       :: ifil
  integer       :: i

! Local variables:

#ifdef MPI
  integer       :: err_ret_code ! value returned by mpi calls 
                                ! (but these mpi calls don't return...)
#endif

  if (ifil .ne. 0) then
    close(unit = ifil)
  end if

#ifdef MPI
! i .gt. 0 implies an error condition, therefore we want to
! kill all the nodes.  This is accomplished with mpi_abort.  If
! it is not an error, exit gracefully with the mpi_finalize.
! NOTE: no mpi functions may be called after a call to mpi_finalize.

  call mpi_group_free(world_group, err_ret_code)
  if (i .eq. 0) then
    call mpi_finalize(err_ret_code)
  else
    call mpi_abort(mpi_comm_world, i, err_ret_code)
  end if
#endif /* MPI */

  if (i .eq. 0) then
    stop
  else
    stop 'PMEMD Terminated Abnormally!'
  end if

end subroutine mexit

!*******************************************************************************
!
! Subroutine:  alloc_error
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_error(routine, string1)

  implicit none
  
  character(*) routine, string1
    
  write(mdout, '(1x,2a)') &
    'FATAL dynamic memory allocation error in subroutine ', routine
  
  write(mdout, '(1x,a)') string1
  
  call mexit(6, 1)
  
end subroutine alloc_error

!*******************************************************************************
!
! Subroutine:  
!
! Description: setup_alloc_error
!              
!*******************************************************************************

subroutine setup_alloc_error

  implicit none
  
  write(mdout, '(1x,2a)') 'FATAL global dynamic memory setup allocation error!'
 
  call mexit(6, 1)

end subroutine setup_alloc_error

end module pmemd_lib_mod
