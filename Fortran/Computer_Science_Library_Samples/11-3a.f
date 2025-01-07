C**********************************************************************
C 11.3�@���d�ϕ��i��`�����j
C**********************************************************************
      implicit none
      integer :: i,j,m,n
      real :: a,b,c,d,h,k,S,f,PAI

      PAI=4.0*atan(1.0)
      m=40  !x������m����
      n=30  !y������n����
      a=-2.0  !x�����̐ϕ����a�`b
      b=2.0
      c=0.0   !y�����̐ϕ����c�`d
      d=PAI
      h=(b-a)/real(m) !x�������ݕ�
      k=(d-c)/real(n) !y�������ݕ�

      S=0.0

      !�ϕ��̈��
      do j=1,n-1
        do i=1,m-1
          S=S+4.0*f(a+i*h,c+j*k)
        end do
      end do

      !4���_���������E
      do i=1,m-1
        S=S+2.0*f(a+i*h,c)
      end do
      do i=1,m-1
        S=S+2.0*f(a+i*h,d)
      end do
      do j=1,n-1
        S=S+2.0*f(a,c+j*k)
      end do
      do j=1,n-1
        S=S+2.0*f(b,c+j*k)
      end do

      !���E��4���_
      S=S+1.0*f(a,c)
      S=S+1.0*f(a,d)
      S=S+1.0*f(b,c)
      S=S+1.0*f(b,d)

      S=h*k/4.0*S
      write(*,901) "S(trapezoid formula) =        ", S

C**********************************************************************
C [�Q�l] ��ϕ��ɂ��v�Z
      write(*,*) " "
      write(*,*) "[for reference]"
      S=88.0/3.0*PAI
      write(*,901) "S(integral from a-b, c-d) =   ", S

  901 format(A30,f15.6)
      end

C**********************************************************************
C �֐�f�̒�`
C**********************************************************************
      real function f(x,y)
      implicit none
      real :: x,y

      f=(x+1)**2-cos(3*y)+5.0
      end function
