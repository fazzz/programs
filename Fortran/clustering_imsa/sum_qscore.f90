
program sum_qscore
  !  use qsort

  implicit none

  integer :: i, j
  integer :: Nmsa ! number of MSAs
  double precision :: d_x, d_y, d_xy

  double precision, allocatable, dimension(:) :: Score_x, Score_y, Score_xy
  
  double precision, allocatable, dimension(:,:) :: d_matrix ! (Nmsa x Nmsa), distance matrix of MSAs

  character inputfilename*200
  
  character(30) :: argv, programname
  character*100 :: arg
  integer :: argc, iargc

  argc=iargc()
  call getarg(0,programname)
  if (argc < 2) then
     call usage(programname)
  else
     call getarg(1, arg)
     read(arg,*), Nmsa
     call getarg(2, inputfilename)
  end if

  open(21,file=inputfilename,status='old')
  allocate(d_matrix(Nmsa,Nmsa))
  call read_dat(21, Nmsa, d_matrix)
  close(21)

  allocate(Score_x(Nmsa), Score_y(Nmsa), Score_xy(Nmsa))
  Score_x = 0.0d0
  Score_y = 0.0d0
  Score_xy = 0.0d0
  
  do i = 1,Nmsa,1
     do j = 1,Nmsa,1
        d_x = d_matrix(i,j)
        d_y = d_matrix(j,i)
        d_xy = 0.5 * ( d_x + d_y )

        Score_x(i) = Score_x(i) + d_x
        Score_y(i) = Score_y(i) + d_y
        Score_xy(i) = Score_xy(i) + d_xy
     end do
  end do

  do i = 1,Nmsa,1
     write(*,'(I3,1X,F7.3,1X,F7.3,1X,F7.3)'), i-1, Score_x(i), Score_y(i), Score_xy(i)
  end do

  deallocate(Score_x,Score_y,Score_xy)

contains
  
  ! Show help mesage
  subroutine usage(programname)
    implicit none
    character(20), intent(in) :: programname
    print *,"Usage: ",programname," Nmsa(Number of MSAs) inputfilename(distance matrix)"
    
    stop
  end subroutine usage
  
  ! Read input file
  subroutine read_dat (unitnum, ndim, x)

    implicit none
    integer, intent(in) :: unitnum
    integer, intent(in) :: ndim
    double precision,allocatable, dimension(:,:), intent(inout) :: x

    integer :: i, j, k
    double precision :: v
    double precision, allocatable, dimension(:) :: A, B
    integer :: err

!    allocate(A(1))
    
    i = 0
    do
       read (unitnum, '(F5.3,1X)', advance="NO",iostat=err, end=999), v

       if (err .eq. 0 ) then
          if ( i > 0) then
             allocate(B(i))
             B = A
             deallocate(A)
          end if
          i = i + 1
          allocate(A(i))
          if ( i - 1 > 0) then
             A(1:i-1) = B(1:i-1)
             deallocate(B)
          end if
          A(i) = v
!          write (*,*) A(i)

       end if
    end do
999 k = 0

    do i=1,ndim,1
       do j=1,ndim,1
          k = k + 1
          x(i,j) = A(k)
       end do
    end do

    deallocate(A)
    
  end subroutine read_dat
  
end program sum_qscore
