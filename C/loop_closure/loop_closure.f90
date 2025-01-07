!===========================================================================
! main program for loop closure 
!===========================================================================
module constants 
  !
  implicit none
  real(8) :: q(3), q_length, q_norm
  real(8) :: qx, qy, qz(2)
  real(8) :: T_alpha(3,3), T_beta(3,3), e1(3)
  real(8),parameter :: PI=3.141592
  real(8),parameter :: alpha=(9.0/180.)*PI,beta=(70.5/180.)*PI
  real(8) :: phi1,r(3),w,test1(3),test2(3),test3(3,3),test4(3,3),theta


  real(8) :: sin1, cos1, sin4(2), cos4(2)
  real(8) :: R1(3,3),R2(3,3,2),R3(3,3,2,2),R4(3,3,2)
  real(8) :: f1,f2,f3,f4
!
contains
!===========================================================================
  subroutine initialize 
    q(1)=3.519 ; q(2)=1.436 ; q(3)=0.
    e1(1)=1.   ; e1(2)=0.   ; e1(3)=0.
    q_norm = dot_product(q, q)
    q_length = dsqrt(q_norm)
    
    T_alpha(1,1)=dcos(alpha) ;T_alpha(1,2)=-dsin(alpha) ;T_alpha(1,3)=0.
    T_alpha(2,1)=dsin(alpha) ;T_alpha(2,2)=dcos(alpha)  ;T_alpha(2,3)=0.
    T_alpha(3,1)=0.          ;T_alpha(3,2)=0.           ;T_alpha(3,3)=1.
    
    T_beta(1,1)=dcos(beta)   ;T_beta(1,2)=-dsin(beta)   ;T_beta(1,3)=0.
    T_beta(2,1)=dsin(beta)   ;T_beta(2,2)=dcos(beta)    ;T_beta(2,3)=0.
    T_beta(3,1)=0.           ;T_beta(3,2)=0.            ;T_beta(3,3)=1.
  end subroutine initialize
!===========================================================================
  subroutine rot_mat_odd(A,theta)   !--- odd => even coordinate rotation
    real(8),intent(inout) :: A(3,3)
    real(8),intent(in) :: theta
    A(1,1)=1.  ; A(1,2)=0.            ; A(1,3)=0.
    A(2,1)=0.  ; A(2,2)=dcos(theta)  ; A(2,3)=-dsin(theta)
    A(3,1)=0.  ; A(3,2)=dsin(theta)  ; A(3,3)=dcos(theta)
  end subroutine rot_mat_odd
!===========================================================================
  subroutine rot_mat_even(A,theta)  !--- even => odd coordinate rotaton
    real(8),intent(inout) :: A(3,3)
    real(8),intent(in) :: theta 
    A(1,1)=1.  ; A(1,2)=0.           ; A(1,3)=0.
    A(2,1)=0.  ; A(2,2)=dcos(theta)  ; A(2,3)=-dsin(theta)
    A(3,1)=0.  ; A(3,2)=dsin(theta)  ; A(3,3)=dcos(theta)
  end subroutine rot_mat_even
!==========================================================================
  subroutine calc_R1(phi1_rad)
    real(8), intent(in) :: phi1_rad 
    call rot_mat_odd(R1,phi1_rad)
  end subroutine calc_R1
!==========================================================================
    subroutine calc_R2(phi1_rad, s, err)
    real(8), intent(in) :: phi1_rad
    real(8), intent(in) :: s(3)
    logical, intent(out) :: err
    integer :: i,j 
    real(8) :: rsquare
    real(8) :: a, D
    real(8) :: sin2, cos2, theta 

    ! -- eq. 21
        r = matmul(transpose(T_beta),matmul(transpose(R1),matmul(transpose(T_alpha),s - q)))
    !=============================for-debug==========================
    !  print *, dsqrt( sum(r(:)**2) ), dsqrt( sum((s(:)-q(:))**2 ) )
    !================================================================

    rsquare = dot_product(r, r)

    ! -- eq. 30
    w = (rsquare - 2.*q(1)*r(1)) / (2. * q(2))
    
    ! -- eq. 29
    a = rsquare - r(1)**2
    D = a - w ** 2

    ! -- eq. 46 
!    if(r(1)/q_length<=f2.or.f1<=r(1)/q_length) then 
!      err = .true.
    if(D<0.) then
      err = .true.
      !=========for-debug===========
      !        print  *, phi1,D
      !=============================
     return
    endif

    ! eq. 23 in paper
    qx=r(1)-q(1)
    qy=w-q(2)
    qz(1)=dsqrt(r(2)**2+r(3)**2-w**2)
    qz(2)=-qz(1)
    !=================for-debug=================
    ! print *, dsqrt(qx**2 + qy**2 + qz(1)**2 )
    !===========================================

    cos2 = ( r(2) * w +  r(3)*dsqrt(D) )/ a
    sin2 = (r(3)*w - r(2)*dsqrt(D))/ a
    theta=datan2(sin2,cos2)
    call rot_mat_even(R2(:,:,1), theta)

    !==========for-debug==============
    !    print *, sin2**2+cos2**2
    !=================================   

    cos2 = ( r(2) * w - r(3)*dsqrt(D) )/ a
    sin2 = (r(3)*w + r(2)*dsqrt(D))/ a
    theta=datan2(sin2,cos2)
    call rot_mat_even(R2(:,:,2),theta )

    !=========for-debug==============
    !    print *, sin2**2+cos2**2
    !================================

    err = .false.
  end subroutine calc_R2
!=============================================================================
  subroutine calc_R4(err)
    logical, intent(out) :: err
    real(8) :: condition,theta(2)
    integer :: i

    ! -- eq. 33
    condition = (q(1) * dcos(beta) &
         - (r(1) - q(1)) * dcos(alpha) &
        - (w    - q(2)) * dsin(alpha) ) &
         / (q(2) * dsin(beta))
    ! -- condition eq.50 
!    if(r(1)/q_length<=f4.or.f3<=r(1)/q_length) then 
!       err = .true.
!       return
!    end if
    if(abs(condition)>=1.) then
       err = .true.
       !====for-debug========
       !  print *, condition
       !=====================
       return
    endif

    ! -- eq. 33
    cos4(1) = condition
    cos4(2) = condition

    ! -- eq. 35
    sin4(1) = dsqrt(1 - condition ** 2)
    sin4(2) = - sin4(1)
    !============for-debug==============
    !     print *, sin4**2+cos4**2
    !===================================

    do i=1,2
       theta(i)=datan2(sin4(i),cos4(i))
       call rot_mat_even(R4(:,:,i), theta(i))
    end do

    err = .false.
  end subroutine calc_R4
!=============================================================================
  subroutine calc_R3
    integer :: i, j
    real(8) :: denominator
    real(8) :: sin3, cos3, theta


    do i = 1, 2
       denominator = (- qx * dsin(alpha) + qy * dcos(alpha)) ** 2 + qz(i) ** 2
       do j = 1, 2
          ! -- eq. 36
          ! -- FIXME: combine intermediate
          cos3 = &
               ((q(1) * dsin(beta) + q(2) * cos4(j) * dcos(beta)) &
                 * (- qx * dsin(alpha) + qy * dcos(alpha)) &
                 + q(2) * qz(i) * sin4(j)) &
                 / denominator
          sin3 = ((q(1) * dsin(beta) + q(2) * cos4(j) * dcos(beta)) * qz(i) &
               - ( -qx * dsin(alpha) + qy * dcos(alpha)) * q(2) * sin4(j)) &
               / denominator

        !============for-debug================
        !    print *, sin3**2+cos3**2                 
        !=====================================

          theta=datan2(sin3,cos3)
          call rot_mat_odd(R3(:, :, i, j), theta)
       end do
    end do
  end subroutine calc_R3
!=============================================================================
    function g( u, R1, R2, R3, R4 ) result(res)
    real(8) :: res
    real(8), intent(in) :: u(3)
    real(8), intent(in) :: R1(3,3), R2(3,3), R3(3,3), R4(3,3)
    real(8) :: prod1(3,3), prod2(3,3), prod3(3,3), prod4(3,3)
    real(8) :: prod_final(3,3), dumy_vec(3),theta , dumy_mat(3,3)

    ! -- eq. 20
    prod1 = matmul(R1, T_beta)
    theta=datan2(R2(3,2),R2(2,2))
    call rot_mat_even(dumy_mat,theta)
    prod2 = matmul(dumy_mat, T_alpha)
    theta=datan2(R4(3,2),R4(2,2)) 
    call rot_mat_even(dumy_mat,theta)
    prod3 = matmul(R3, T_beta)
    prod4 = matmul(dumy_mat, T_alpha)

    prod_final = &
         matmul(T_alpha, &
         matmul(prod1, &
         matmul(prod2, &
         matmul(prod3, &
         prod4))))           

   dumy_vec = matmul(prod_final, e1)
   res = dot_product(u, dumy_vec) - dcos(beta)

  end function g
!==========================================================================
  subroutine boundary(d)
    real(8) :: d, dsquare,dumy1,dumy2
  
    ! -- eq. 47 in paper
    dsquare = d ** 2

    dumy1 = q(1) * dsquare
    dumy2 = q(2) * d * dsqrt(4.*q_norm - dsquare)
    f1 = (dumy1+dumy2) / (2.*q_norm*q_length)
    f2 = (dumy1-dumy2) / (2.*q_norm*q_length)

    ! -- eq. 51
    dumy1 = 2.*q(2)*q(1)*(dcos(alpha)+dcos(beta)) - dsquare*dsin(alpha)
    dumy2 = 2. * (q(2)*dcos(alpha) - q(1)*dsin(alpha))

    f3=(2*(q(2)**2)*(dsin(alpha)+dsin(beta)) + dumy1)/(dumy2*q_length)
    f4=(2*(q(2)**2)*(dsin(alpha)-dsin(beta)) + dumy1)/(dumy2*q_length)

!===========================for debug============================
!    write(10,'(e16.10,4(1x,e16.10))') (d/q_length),f1,f2,f3,f4
!================================================================
  end subroutine boundary
!========================================================================
end module constants
!==========================================================================

!==========================================================================
!==========================================================================
program main
  use constants
  integer :: i,j,Num
  real(8) :: s(3),u(3)
  real(8) :: range_r
  real(8) :: rsquare
  real(8),allocatable :: P(:,:)
  real(8) :: rr,r_length
  real(8) :: phi1_rad
  logical :: err
  !
  open(30,file='f.dat')
  !
  call initialize

! =============for debug ==============
!  open(10,file='boundary.dat')
!  do range_r=0. , 2.*q_length , 0.001
!     call boundary(range_r) 
!  enddo
!  stop
!======================================
 
  open(20,file='C_alpha_coor_alpha_helix_resi5.dat')
  read(20,*) Num
  allocate(P(Num,3)) 
  read(20,*) u(3)
!======= for debug========
!   u(1)=0.730;u(2)=-0.280;u(3)=0.624
!  u(1)=0.;u(2)=1.;u(3)=0.
!  u(1)=0.059;u(2)=-0.852;u(3)=-0.521
!=========================
  do i=1,Num
     read(20,*) P(i,:)
  enddo
  s(:)=P(Num-1,:)-P(1,:)
!====== for debug=========
!  s(1)=2.913;s(2)=-3.398;s(3)=-2.407
! s(1)=9.;s(2)=0.;s(3)=0.
! s(1)=8.821;s(2)=1.824;s(3)=-3.101
!=========================

  r_length=dsqrt(sum((s(:)-q(:))**2))
  temp1=r_length/q_length


!  if( temp1<=1.2051.or.1.9069<=temp1 ) then
!    print *,  'r/q =' , temp1 , 'ERROR! r/q should be 1.2051<=r/q<=1.9069'
!     stop
!  endif

  call boundary(r_length)
!===== for debug =================
!  print *, 'r/q_length=' , temp1
!  print *, f2, '<=x/q<=' , f1
!  print *, f4, '<=x/q<=' , f3
!=================================

  do phi1 = 0. , 360.,0.001
     phi1_rad = (phi1/180.) * PI 

     call calc_R1(phi1_rad)
     call calc_R2(phi1_rad, s, err)

     if (err .eqv. .true.) then
       cycle
     endif

     ! -- R4 should be calculated before R3
     call calc_R4(err)
     if (err .eqv. .true.) then
       cycle
     endif

     call calc_R3
     
     do i = 1, 2
        do j = 1, 2
           write(30, *) phi1, g( u(:), R1(:,:) , R2(:,:,i), R3(:,:,i,j), R4(:,:,j))
!           test1(:)=matmul( transpose(T_beta) ,matmul( transpose(R1), matmul( transpose(T_alpha) , s - q )))
!            print *, phi1
!            print *, test1(:)
!            print *, matmul(R2(:,:,i),q) + matmul(R2(:,:,i),matmul(T_alpha,matmul(R3(:,:,i,j),matmul(T_beta,matmul(R4(:,:,j),q)))))
!            test2(1)=qx ; test2(2)=qy ; test2(3)=qz(i)
!            print *, matmul(R2(:,:,i),q+test2)
!           print *, matmul(transpose(T_beta),matmul(transpose(R3(:,:,i,j)),matmul(transpose(T_alpha),test2)))
!           print *, q
!           print *, R3(1,1:3,i,j)
!           print *, R3(2,1:3,i,j)
!           print *, R3(3,1:3,i,j)
!           print *, R4(1,1:3,j)
!           print *, R4(2,1:3,j)
!           print *, R4(3,1:3,j)
!           print *, T_beta(1,1:3)
!           print *, T_beta(2,1:3)
!           print *, T_beta(3,1:3)
!           theta=datan2(R2(3,2,i),R2(2,2,i))+PI
!           call rot_mat_even(test3,theta)
!           theta=datan2(R4(3,2,j),R4(2,2,j))+PI
!           call rot_mat_even(test4,theta)
!           print *, matmul(R1(:,:),matmul(R2(:,:,i),matmul(T_beta,matmul(R3(:,:,i,j),matmul(R4(:,:,j),q)))))
!           print *, matmul(R1(:,:),matmul(test3,matmul(T_beta,matmul(R3(:,:,i,j),matmul(test4,q)))))
!           print *, matmul(R2(:,:,i),matmul(T_beta,matmul(R4(:,:,j),q)))
!           print *, matmul(test3,matmul(T_beta,matmul(test4,q)))
!           test1(1)=qx ; test1(2)=qy ; test1(3)=qz(i)
!           test2(:)=matmul( transpose(R3(:,:,i,j)) ,matmul( transpose(T_beta), matmul( transpose(T_alpha) , s - q )))
!           test2(1)=qx ; test2(2)=qy ; test2(3)=qz(i)
!           print *, test1(:)-matmul(R2(:,:,i),q+test2)
!           test(1)=qx ; test(2)=qy ; test(3)=qz(i)
!           print *, phi1,real(matmul(transpose(R3(:,:,i,j)),matmul(transpose(T_alpha),test))-matmul(T_beta,matmul(R4(:,:,j),q)))
!           test2(:)=matmul(T_alpha,matmul(R3(:,:,i,j),matmul(T_beta,matmul(R4(:,:,j),q))))
!           print *, test(:)-matmul(T_alpha,matmul(R3(:,:,i,j),matmul(T_beta,matmul(R4(:,:,j),q))))
!           print *, 'c',matmul(R2(:,:,i),q+test)
!          print *, matmul(T_beta,matmul(R4(:,:,j),q))
        enddo
     enddo
  enddo
end program main
!=========================================================================================
!=========================================================================================
