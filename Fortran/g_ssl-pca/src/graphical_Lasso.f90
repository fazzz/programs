! This module contains functions / subroutines 
! which conduct Graphical Lasso calculatin.

!  subroutine set_W(W,omg,sig,iA,m,i)   : set W matrix( W, omega,sig <- iA )
!  subroutine set_S(matS,vecS,eleS,m,i) : set vector s and s(i,i) from matrix S ( VC matrix)
!  subroutine set_A(A,vecl,lambda,m,i)  : set A from vecl and lambda
!  subroutine set_iA(iA,vecS,W,beta,vecA,m,i)  : set iA from W, beta, and vecA
!  subroutine Lasso_Reg(beta,A,rou,vecS,W,m)   : calc. of beta and A using iterativemethod
!  subroutine block_coordinate_descent_method(A,iA,S,rou,m) : A,iA <- S,rou.m
!  subroutine cal_anomaly_x_bet_A(xn,A,an,m) : calc. anomaly of xn from matrix A

module graphical_Lasso
  implicit none

contains

  subroutine set_W(W,omg,sig,iA,m,i) ! set W matrix( W, omega, sig <- iA )
    implicit none

    double precision, allocatable, dimension(:,:), intent(in)  ::   iA  ! Precision Matrix (m,m)
    double precision, allocatable, dimension(:,:), intent(inout) :: W   ! W matrix (m-1,m-1)
    double precision, allocatable, dimension(:),   intent(inout) :: omg ! omega (m-1)
    double precision, intent(out) :: sig  ! sigma
    integer, intent(in) :: m ! m : # of dimensions 
    integer, intent(in) :: i ! i : specified element number "i"

    integer j,k

    if ( i /= 1 .and. i /= m ) then
       W(1:i-1,1:i-1) = iA(1:i-1,1:i-1)
       W(i:m-1,1:i-1) = iA(i+1:m,1:i-1)
       W(1:i-1,i:m-1) = iA(1:i-1,i+1:m)
       W(i:m-1,i:m-1) = iA(i+1:m,i+1:m)

       omg(1:i-1) = iA(i,1:i-1)
       omg(i:m-1) = iA(i,i+1:m)
    else if ( i == 1 ) then
       W = iA(2:m,2:m)
       omg(1:m-1) = iA(1,2:m)
    else
       W = iA(1:m-1,1:m-1)
       omg(1:m-1) = iA(m,1:m-1)
    endif

    sig = iA(i,i)

    ! Outputs for debug
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write(*,'("W=")')                            !
    ! do j=1,m-1,1                                 !
    !    do k=1,m-1,1                              !
    !       write(*,'(F8.3)',advance='no'), W(j,k) !
    !    end do                                    !
    !    write(*,'(/)')                            !
    ! end do                                       !
    !                                              !
    ! write(*,'("omg=")')                          !
    ! do j=1,m-1,1                                 !
    !    write(*,'(F8.3)',advance='no'),omg(j)     !
    ! end do                                       !
    ! write(*,'(/)')                               !
    !                                              !
    ! write(*,'("sig=")',advance='no')             !
    ! write(*,'(F8.3)'),sig                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Outputs for debug

  end subroutine set_W

  subroutine set_S(matS,vecS,eleS,m,i) ! set vector s and s(i,i) from matrix S (VC matrix)
    implicit none

    double precision, allocatable, dimension(:,:), intent(in) :: matS  ! Variance-Covariance mat (m,m)
    double precision, allocatable, dimension(:), intent(inout)  :: vecS  ! vector s (m-1)
    double precision, intent(out) :: eleS ! matS(i,i)
    integer, intent(in) :: m ! m : # of dimensions 
    integer, intent(in) :: i ! i : specified element "i"

    integer j,k

    ! Outputs for debug
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write(*,'("matS=")')                               !
    ! do j=1,m,1                                         !
    !    do k=1,m,1                                      !
    !       write(*,'(F8.3,1X)',advance='no'), matS(j,k) !
    !    end do                                          !
    !    write(*,'(/)',advance='no')                     !
    ! end do                                             !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Outputs for debug

    if ( i /= 1 .and. i /= m ) then
       vecS(1:i-1) = matS(i,1:i-1)
       vecS(i:m-1) = matS(i,i+1:m)
    else if ( i .eq. 1 ) then
       vecS(1:m-1) = matS(1,2:m)
    else
       vecS(1:m-1) = matS(m,1:m-1)
    endif

    eleS = matS(i,i)

    ! Outputs for debug
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write(*,'("vecS=")')          !
    ! do j=1,m-1,1                  !
    !    write(*,'(F8.3 )'),vecS(j) !
    ! end do                        !
    ! write(*,'(/)')                !
    !                               !
    ! write(*,'("eleS=")')          !
    ! write(*,'(F8.3 )'),eleS       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Outputs for debug

  end subroutine set_S

  subroutine set_A(A,vecl,lambda,m,i)  ! A <- vecl, lambda
    implicit none

    double precision, allocatable, dimension(:), intent(in) :: vecl    ! vector l (m-1)
    double precision, intent(in) ::             lambda                 ! lambda
    double precision, allocatable, dimension(:,:), intent(inout) ::  A ! Precision Matrix (m,m)
    integer, intent(in) :: m ! m : # of dimensions 
    integer, intent(in) :: i ! i : specified element "i"

    integer j,k

    if ( i /= 1 .and. i /= m ) then
       A(i,1:i-1) = vecl(1:i-1)
       A(i,i+1:m) = vecl(i:m-1)

       A(1:i-1,i) = vecl(1:i-1)
       A(i+1:m,i) = vecl(i:m-1)
    else if ( i .eq. 1 ) then
       A(1,2:m) = vecl(i:m-1)
       A(2:m,i) = vecl(i:m-1)
    else
       A(m,1:m-1) = vecl(1:m-1)
       A(1:m-1,m) = vecl(1:m-1)
    endif

    A(i,i) = lambda

    ! Outputs for debug
    !    write(*,'("A=")')
    !    do j=1,m,1
    !       do k=1,m,1
    !          write(*,'(F8.3)',advance='no'),A(j,k)
    !       end do
    !       write(*,'(/)')
    !    end do
    ! Outputs for debug

  end subroutine set_A

  subroutine set_iA(iA,omg,sig,m,i)  ! iA <- vecl, sig
    implicit none

    double precision, allocatable, dimension(:), intent(in) :: omg   ! vector omega (m-1)
    double precision, intent(in) ::                sig               ! sigma
    double precision, allocatable, dimension(:,:), intent(inout) :: iA ! Precision Matrix (m,m)
    integer, intent(in) :: m ! m : # of dimensions 
    integer, intent(in) :: i ! i : specified element "i"

    integer j, k
  
    if ( i /= 1 .and. i /= m ) then
       iA(i,1:i-1) = omg(1:i-1)
       iA(i,i+1:m) = omg(i:m-1)

       iA(1:i-1,i) = omg(1:i-1)
       iA(i+1:m,i) = omg(i:m-1)
    else if ( i .eq. 1 ) then
       iA(1,2:m) = omg(1:m-1)
       iA(2:m,1) = omg(1:m-1)
    else
       iA(m,1:m-1) = omg(1:m-1)
       iA(1:m-1,m) = omg(1:m-1)
    endif

    iA(i,i) = sig

    ! Outputs for debug
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write(*,'("sig=",F8.3)'),sig                 !
    ! write(*,'("i=",I4)'),i                       !
    ! write(*,'("iA=")')                           !
    ! do j=1,m,1                                   !
    !    do k=1,m,1                                !
    !       write(*,'(F8.3)',advance='no'),iA(j,k) !
    !    end do                                    !
    !    write(*,'(/)')                            !
    ! end do                                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Outputs for debug

  end subroutine set_iA

  subroutine Lasso_Reg(beta,A,rou,vecS,W,m) ! calc. of beta and A using iterative method
    implicit none

    double precision, allocatable, dimension(:,:), intent(in) :: W    ! matrix W (m-1,m-1)
    double precision, allocatable, dimension(:), intent(in) ::  vecS  ! vector S (m-1)
    double precision, intent(in) :: rou                               ! the degree of sparsity
    double precision, allocatable, dimension(:), intent(inout) ::   A   ! vector A (m-1)
    double precision, allocatable, dimension(:), intent(inout) ::  beta ! vector beta (m-1)
    double precision, allocatable, dimension(:) ::  beta_prv          ! vector beta_prv (m-1)
    double precision, allocatable, dimension(:) ::  beta_dif          ! vector beta_dif (m-1)
    double precision  max_diff, min_diff                              ! value for convergence check
    integer, intent(in) :: m ! m : # of dimensions

    double precision, parameter :: torelance = 1.0e-5

!    integer i,j,k
    integer j,k
    integer num_ite

    num_ite = 1

    A(1:m-1) = vecS(1:m-1)

    max_diff = 0.0d0

    do while (max_diff > torelance .or. num_ite < 3)  ! converegence of beta
       do j=1,m-1,1
          if (A(j) > rou)  then
             beta(j) = (A(j) - rou) / W(j,j)
          else if (A(j) < -rou) then
             beta(j) = (A(j) + rou) / W(j,j)
          else
             beta(j) = 0.0d0
          endif
       end do

       do j=1,m-1,1
          A(j) = vecS(j)
          do k=1,m-1,1
             if (k /= j) then
                A(j) = A(j) - W(j,k) * beta(k)
             endif
          end do
       end do

       if ( num_ite > 1) then
          beta_dif = beta - beta_prv
          max_diff = maxval(beta_dif)
          max_diff = abs(max_diff)
          min_diff = minval(beta_dif)
          min_diff = abs(min_diff)
          if ( min_diff > max_diff ) then
             max_diff = min_diff
          end if
       endif

       beta_prv = beta
!       if ( num_ite > 1) then
!          write(*,'(I5,"-th max diff:",F10.6)'),num_ite,max_diff
!       end if
       num_ite = num_ite + 1

    end do

    ! Outputs for debug
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write(*,'("beta=")')         !
    ! do j=1,m-1,1                 !
    !    write(*,'(F8.3)'),beta(j) !
    ! end do                       !
    ! write(*,'("A=")')            !
    ! do j=1,m-1,1                 !
    !    write(*,'(F8.3)'),A(j)    !
    ! end do                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Outputs for debug

  end subroutine Lasso_Reg

  subroutine block_coordinate_descent_method(A,iA,S,rou,m) ! A,iA <- S,rou.m
    implicit none

    double precision, allocatable, dimension(:,:), intent(in) ::       S  ! Variance-Covariance Mat. (m,m)
    double precision, allocatable, dimension(:,:), intent(inout) :: A, iA ! Precision Matrix (m,m)
    double precision, intent(in) :: rou                                   ! degree of the sparsity
    integer, intent(in) :: m ! m : # of dimensions 

    double precision, allocatable, dimension(:,:) ::   W   ! W matrix (m-1,m-1)
    double precision, allocatable, dimension(:)   ::   omg ! omega (m-1)
    double precision  sig  ! sigma

    double precision, allocatable, dimension(:)   :: vecS  ! vector s (m-1)
    double precision, allocatable, dimension(:)   :: vecA  ! vector A (m-1)
    double precision, allocatable, dimension(:)   :: vecl  ! vector l (m-1)
    double precision eleS ! matS(i,i)
    double precision lambda

    double precision, allocatable, dimension(:)   ::  beta ! vector beta (m-1)

    double precision, allocatable, dimension(:,:) :: AiA
    double precision, allocatable, dimension(:)   :: Wbeta      ! matrix W x vector beta (m-1)
    double precision                              :: betaTWbeta ! vector betaT matrix W x vector beta

    double precision, allocatable, dimension(:,:) :: A_prv, diff_A 
    double precision  max_diff, min_diff                        ! values for convergence check
    double precision, parameter :: torelance = 1.0e-5

    integer i,j,k
    integer num_ite

    allocate(W(m-1,m-1))
    allocate(omg(m-1))
    allocate(vecS(m-1))

    allocate(vecA(m-1))
    allocate(beta(m-1))
    allocate(Wbeta(m-1))

    allocate(vecl(m-1))

    allocate(A_prv(m,m))

    allocate(AiA(m,m))

    !    write(*,'("block_coordinate_descent_method")')

    iA = S    ! initialization
    A = 0.0d0 ! initialization

    num_ite = 1
    max_diff = 0.0e0

    do while(max_diff > torelance .or. num_ite < 3) ! Converegence of A and A^-1
       A_prv = iA
       do i=1,m,1
          call set_S(S,vecS,eleS,m,i)      ! vecS, eleS  <- matS, i, m
          call set_W(W,omg,sig,iA,m,i)        ! W, omg, sig <- iA, i, m

          call Lasso_Reg(beta,vecA,rou,vecS,W,m) ! beta, vecA, <- rou, vecS, W, m

          ! omega = W * beta
          omg = 0.0d0
          do j=1,m-1,1
             do k=1,m-1,1
                omg(j) = omg(j) + W(j,k) * beta(k)
             end do
          end do

          sig = eleS + rou

          ! Wbeta = W * beta
          Wbeta = 0.0d0
          do j=1,m-1,1
             do k=1,m-1,1
                Wbeta(j) = Wbeta(j) + W(j,k) * beta(k)
             end do
          end do

          ! betaTWbeta = betaT * Wbeta
          betaTWbeta = 0.0d0
          do j=1,m-1,1
             betaTWbeta = betaTWbeta + beta(j) * Wbeta(j)
          end do

          lambda = 1.0d0 / (sig - betaTWbeta)
          do j=1,m-1,1
             vecl(j) = -(beta(j)*lambda)
          end do

          call set_iA(iA,omg,sig,m,i)    ! iA <-  omg,    sig, m

          call set_A(A,vecl,lambda,m,i)  !  A <- vecl, lambda, m

       end do

       if ( num_ite > 1) then
          diff_A = iA - A_prv
          max_diff = maxval(diff_A)
          max_diff = abs(max_diff)
          min_diff = minval(diff_A)
          min_diff = abs(min_diff)
          if ( min_diff > max_diff ) then
             max_diff = min_diff
          end if
!          write(*,'(I3,"-th max diff:",F10.6)'),num_ite,max_diff
       else
!          write(*,'("The begining of iteration of block coordinate descent method... ")')
       endif

       num_ite = num_ite + 1
    end do

    ! Outputs for debug
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write (*,'("A=")')                               !
    ! do i=1,m                                         !
    !    do j=1,m                                      !
    !       write (*,'(F8.3,1X)',advance='no'),A(i,j)  !
    !    end do                                        !
    !    write (*,'(/)')                               !
    ! end do                                           !
    ! write (*,'("iA=")')                              !
    ! do i=1,m                                         !
    !    do j=1,m                                      !
    !       write (*,'(F8.3,1X)',advance='no'),iA(i,j) !
    !    end do                                        !
    !    write (*,'(/)')                               !
    ! end do                                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Outputs for debug

    AiA = 0.0d0
    do i=1,m
       do j=1,m
          do k=1,m
             AiA(i,j) = AiA(i,j) + A(i,k)*iA(k,j)
          end do
       end do
    end do

    ! Outputs for debug
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! write (*,'("A*iA=")')                             !
    ! do i=1,m                                          !
    !    do j=1,m                                       !
    !       write (*,'(F8.3,1X)',advance='no'),AiA(i,j) !
    !    end do                                         !
    !    write (*,'(/)')                                !
    ! end do                                            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Outputs for debug

!    write(*,'("The block coordinate descent method is converged.")')
    
    deallocate(W)
    deallocate(omg)

    deallocate(vecA)
    deallocate(beta)

    deallocate(AiA)

  end subroutine block_coordinate_descent_method

  subroutine cal_anomaly_x_bet_A(xn,A,an,m) ! calc. anomaly of xn from matrix A
    implicit none

    double precision, allocatable, dimension(:,:), intent(in) :: A ! Precision Matrix (m,m)
    double precision, allocatable,dimension(:), intent(in) :: xn    ! (m), m dimension
    double precision, allocatable,dimension(:), intent(inout) :: an    ! (m), m dimension
    integer, intent(in) :: m ! m : # of dimensions 

    integer i,j
    double precision dpi, sum
    double precision, allocatable,dimension(:) :: invA ! (m)

    dpi = acos(-1.0d0)
    
    an = 0.0d0

    allocate(invA(m))
    
    do i = 1,m,1
       invA(i) = 1.0d0 / A(i,i)
    end do
       
    do i = 1,m,1
       an(i) = 0.5d0 * log(dpi) * invA(i)
       sum = 0.0d0
       do j = 1,m,1
          sum = sum + A(i,j)*xn(j)
       end do
       an(i) = an(i) + 0.5 *invA(i) * sum * sum
    end do

    deallocate(invA)
    
  end subroutine cal_anomaly_x_bet_A

end module graphical_Lasso
