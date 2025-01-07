! This module contains the subroutines to perform
! Principle Component Analysis (PCA)

module pca
  implicit none

contains
  ! Calculate center of mass
  subroutine get_com(natom, frame, coord, mass, com)
    implicit none

    integer :: i, j
    integer, intent(in) :: natom, frame
    double precision, intent(in)  :: coord(3, natom, frame), mass(natom)
    double precision, intent(out) :: com(3, frame)
    double precision :: summass
    
    summass = sum(mass)
    if(summass == 0) write(*,*), "error: sum of mass = 0"
    do j = 1, frame
       do i = 1, 3
          com(i, j) = dot_product(coord(i,:,j), mass(:)) / summass
       end do
    enddo
  end subroutine get_com

  ! Calcurate mean structure
  subroutine bftrj(natom, frame, coord, mass, com, bfcoord)
    use bestfit

    implicit none

    integer :: i, j, k, l, m
    integer, intent(in) :: natom, frame

    double precision, allocatable, dimension(:,:),   intent(inout)  :: com
    double precision, allocatable, dimension(:),     intent(in)  :: mass
    double precision, allocatable, dimension(:,:,:), intent(inout)  :: coord
    double precision, allocatable, dimension(:,:,:), intent(inout) :: bfcoord

    double precision, allocatable, dimension(:,:) :: sumcoord, meancoord, outcoord, bfsnap
    
    call get_com(natom, frame, coord, mass, com)

    allocate(meancoord(3, natom),  &
              sumcoord(3, natom),  outcoord(3, natom), bfsnap(3, natom))

    do j = 1, frame
       do i = 1, natom
          bfcoord(:,i,j) = coord(:,i,j) - com(:,j)
       enddo
    enddo
    
    do k = 1,5,1
       do i = 1, natom
          sumcoord(1,i) = sum(bfcoord(1,i,:))
          sumcoord(2,i) = sum(bfcoord(2,i,:))
          sumcoord(3,i) = sum(bfcoord(3,i,:))
       enddo
       meancoord(:,:) = sumcoord(:,:) / frame

       do i = 1,frame,1
          bfsnap(:,:) = bfcoord(1:3,1:natom,i)

          call fit(natom, meancoord, bfsnap, mass, outcoord)
          bfcoord(:,:,i) = outcoord(:,:)
       enddo
    enddo

    deallocate(meancoord, sumcoord, outcoord, bfsnap)
    
  end subroutine bftrj

  ! Variance-covariance matrix
  subroutine varcor(natom, bfcoord, varm, frame, mass)
    implicit none

    integer :: i, j, k, l, natom, frame
    double precision, allocatable, dimension(:,:,:), intent(in) :: bfcoord
    double precision, allocatable, dimension(:),     intent(in) :: mass
    double precision, allocatable, dimension(:,:),   intent(out) :: varm
    double precision, allocatable, dimension(:,:) :: x
    double precision, allocatable, dimension(:)   :: xm,  nstmass

    allocate(x(natom*3, frame), varm(natom*3, natom*3))
    allocate(xm(natom*3), nstmass(natom*3))

    ! Convert 3*N*frame matrix to 3N*frame matrix
    do i = 1, frame
       x(1:natom,i)           = bfcoord(1,1:natom,i)
       x(natom+1:natom*2,i)   = bfcoord(2,1:natom,i)
       x(natom*2+1:natom*3,i) = bfcoord(3,1:natom,i)
    enddo
    
    nstmass(1:natom)=mass(1:natom)
    nstmass(natom+1:natom*2)=mass(1:natom)
    nstmass(natom*2+1:natom*3)=mass(1:natom)

    ! Average
    do i = 1, natom*3
       xm(i) = 0.0D0
       do j = 1, frame
          xm(i) = xm(i) + x(i,j)*sqrt(nstmass(i))
       enddo
       xm(i) = xm(i) / dble(frame)
    enddo

    ! Variance-covariance
    do j = 1, natom*3
       do k = 1, natom*3
          varm(j,k) = 0.0D0
          do l = 1, frame
             varm(j,k) = varm(j,k)+(x(j,l)-xm(j))*(x(k,l)-xm(k))
          enddo
          varm(j,k) = varm(j,k) / dble(frame)
       enddo
    enddo

    deallocate(x)
    deallocate(xm, nstmass)
    
  end subroutine varcor

  ! Eigen value analysis
  subroutine eigvec(natom, varm, eigv)
    implicit none

    integer :: i, j, k
    integer, intent(in) :: natom
    double precision, allocatable, dimension(:,:), intent(in)  :: varm
    double precision, allocatable, dimension(:), intent(inout) :: eigv

    integer :: n, lda, lwork, info
    double precision, allocatable, dimension(:) :: workspace

    external dsyev
    
    n = natom*3
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
    
  end subroutine eigvec
  
end module pca
