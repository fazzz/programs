! This module standradlize the data, time series .
! x <- x-m/(sqrt(s))

!  subroutine ave_timeseries(x,avx,n,m)     : time average of each dimension.
!  subroutine vcv_timeseries(x,n,m,avx,vcx) : variance covariance matrix
!  subroutine standradlize_ts(x,n,m)        : standradlize the data , time series

module standardlize
  implicit none

contains
  
  subroutine ave_timeseries(x,avx,n,m) ! time average of each dimension.
    implicit none

    double precision, allocatable, dimension(:,:), intent(in)  :: x   ! data (m,n)
    double precision, allocatable, dimension(:), intent(inout) :: avx ! average of data (m)
    integer, intent(in) :: m ! m : # of dimensions 
    integer, intent(in) :: n ! n : # of steps

    integer i,j
    double precision sum

    do i=1,m,1
       sum = 0.0e0
       do j=1,n,1
          sum = sum + x(i,j)
       end do
       avx(i) = sum / n
    end do

  end subroutine ave_timeseries

  subroutine vcv_timeseries(x,n,m,avx,vcx) ! variance covariance matrix
    implicit none

    double precision, allocatable,dimension(:,:),intent(in) :: x  ! data (m,n)
    double precision, allocatable,dimension(:),intent(in) :: avx  ! average of data (m)
    integer, intent(in) :: m ! m : # of dimensions 
    integer, intent(in) :: n ! n : # of steps

    double precision, allocatable,dimension(:,:),intent(inout) :: vcx ! variance of data (m,m)

    integer i,j,k

    do i=1,m,1
       do j=1,m,1
          vcx(j,i) = 0.0e0
          do k=1,n,1
             vcx(j,i) = vcx(j,i) + ( (x(i,k) - avx(i)) * (x(j,k) - avx(j)) )
          end do
          vcx(j,i) = vcx(j,i) / n
       end do
    end do

  end subroutine vcv_timeseries

  subroutine standradlize_ts(x,n,m) ! standradlize the data , time series
    implicit none

    double precision, allocatable,dimension(:,:),intent(inout) :: x ! data(m,n)
    integer,intent(in) :: m ! # of dimensions
    integer,intent(in) :: n ! # of steps

    double precision, allocatable,dimension(:) :: avx ! average of data (m)
    double precision, allocatable,dimension(:,:) :: vcx ! variance-covariance of data (m,m)

    integer i,j

    allocate(avx(m))
    allocate(vcx(m,m))

    call ave_timeseries(x,avx,n,m)

    call vcv_timeseries(x,n,m,avx,vcx)

    do i=1,n,1
       do j=1,m,1
          x(j,i) = (x(j,i) - avx(j)) / sqrt(vcx(j,j))
       end do
    end do

    deallocate(avx)
    deallocate(vcx)

  end subroutine standradlize_ts

end module standardlize


