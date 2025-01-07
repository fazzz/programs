! This modeule contains the function to perform
! Potential Energy Principle Component Analysis-sincos (PEPCA)

module pepca
  implicit none

contains

  ! Variance-covariance matrix
  subroutine enevc(nene, ene, varm, frame)
    implicit none

    integer :: i, j, k, l 
    integer , intent(in) :: nene, frame
    double precision, allocatable, dimension(:,:), intent(in)  :: ene
    double precision, allocatable, dimension(:,:), intent(out) :: varm
    double precision, allocatable, dimension(:,:) :: x
    double precision, allocatable, dimension(:) :: enem

    allocate(varm(nene,nene))
    allocate(enem(nene))
        
    ! Average
    do i = 1, nene
       enem(i) = 0.0D0
       do j = 1, frame
          enem(i) = enem(i) + ene(i,j)
       enddo
       enem(i) = enem(i) / dble(frame)
    enddo

    ! Variance-covariance
    do j = 1, nene
       do k = 1, nene
          varm(j,k) = 0.0D0
          do l = 1, frame
             varm(j,k) = varm(j,k)+(ene(j,l)-enem(j))*(ene(k,l)-enem(k))
          enddo
          varm(j,k) = varm(j,k) / dble(frame)
       enddo
    enddo

    deallocate(enem)
    
  end subroutine enevc

  ! Eigen value analysis
  subroutine pepcaeigvec(nene, varm, eigv)
    implicit none

    integer :: i, j, k
    integer, intent(in) :: nene
    double precision, allocatable, dimension(:,:), intent(inout)  :: varm
    double precision, allocatable, dimension(:), intent(inout) :: eigv

    integer :: n, lda, lwork, info
    double precision, allocatable, dimension(:) :: workspace

    external dsyev
    
    n = nene
    lda=n
    lwork=n*300

    allocate(workspace(lwork))
    
    call dsyev('V', 'u', n, varm, lda, eigv,workspace, lwork, info)

    if (info/=0) Then
       write (6, *) 'Failure in DSYEV.  INFO = ', info
       stop
    endif

    if (lwork < workspace(1)) Then
       write (6, *) 'Optimum workspace required = ', workspace(1), &
            'Workspace provided         = ', lwork
       stop
    endif

    deallocate(Workspace)
    
  end subroutine pepcaeigvec
  
end module pepca
