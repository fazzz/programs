! This program performs simple clustering for multiple MSAs
! Input : Similarity(distance^-1) matrix of MSAs, threshold value for distance
! Output : Indices and Num. of Neiboring MSAs for each MSA

program clust_simple
!  use qsort
  
  implicit none

  integer :: i, j
  integer :: Nmsa ! number of MSAs
  double precision :: d_x, d_y, d_xy
  
  integer, allocatable, dimension(:) :: Numc ! (Nmsa), Number of MSAs which is neibor to a specific MSA
  double precision, allocatable, dimension(:) :: Score
  
  double precision :: d_threshold ! the threshold value for the distance between MSAs
  double precision, allocatable, dimension(:,:) :: d_matrix ! (Nmsa x Nmsa), distance matrix of MSAs

  character inputfilename*200
  
  character(30) :: argv, programname
  character*100 :: arg
  integer :: argc, iargc

  argc=iargc()
  call getarg(0,programname)
  if (argc < 3) then
     call usage(programname)
  else
     call getarg(1, arg)
     read(arg,*), Nmsa
     call getarg(2, arg)
     read(arg,*), d_threshold
     call getarg(3, inputfilename)
  end if

!  write(*,'(" Nmsa:", I4)'), Nmsa
!  write(*,'(" d_threshold:", F4.3)'), d_threshold
  
  open(21,file=inputfilename,status='old')
  allocate(d_matrix(Nmsa,Nmsa))
  call read_dat(21, Nmsa, d_matrix)
  close(21)

!  do i=1,Nmsa,1
!     do j=1,Nmsa,1
!        write (*, '(F5.3,1X)', advance='no'), d_matrix(i,j)
!     end do
!     write (*, '(1X)')
!  end do
  
  allocate(Numc(Nmsa))
  Numc = 0

  allocate(Score(Nmsa))
  Score = 0.0d0
  
  do i = 1,Nmsa,1
!     do j = i+1,Nmsa,1
     do j = 1,Nmsa,1
        d_x = d_matrix(i,j)
        d_y = d_matrix(j,i)
        d_xy = 0.5*( d_x + d_y )
        

        if ( d_xy .gt. d_threshold) then
           Numc(i) = Numc(i) + 1
           Score(i) = Score(i) + d_xy
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! d = d_matrix(j,i)             !
        ! if ( d .gt. d_threshold) then !
        !    Numc(j) = Numc(j) + 1      !
        !    Score(j) = Score(j) + d    !
        ! end if                        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
     end do
  end do

  do i = 1,Nmsa,1
!     write(*,'(I3,1X,I3,1X,F7.3)'), i-1, Numc(i), Score(i)
     write(*,'(F7.3)'), Score(i)
  end do

  deallocate(Numc)
  deallocate(Score)

contains

  ! Show help mesage
  subroutine usage(programname)
    implicit none
    character(20), intent(in) :: programname
    print *,"Usage: ",programname," Nmsa(Number of MSAs) d_threshold(threshold value) inputfilename(distance matrix)"
    
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

end program clust_simple
