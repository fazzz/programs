C**********************************************************************
C 10.1�@�敪���ϖ@�Ƒ�`����
C**********************************************************************
      implicit none
      integer :: n,j
      real :: a,b,f,h,S,xj,xj_1

C �A���S���Y���菇�P�D
      n=10
      a=0.0
      b=2.0*atan(1.0) !��/2

C �A���S���Y���菇�Q�D
      h=(b-a)/real(n)
      S=0.0

C �A���S���Y���菇�R�D
      do j=1,n-1
        S=S+f(a+j*h)
      end do

C �A���S���Y���菇�S�D
      S=(f(a)+2.0*S+f(b))*h/2.0

C ���ʂ̏o�́i��`�����̌��ʁj
      write(*,901) "S(trapezoid formula) =        ", S

C**********************************************************************
C [�Q�l] ��ϕ��ɂ��v�Z
      write(*,*) " "
      write(*,*) "[for reference]"
      S=-cos(b)+cos(a)
      write(*,901) "S(integral from a to b) =     ", S

C [�Q�l] �敪���ϖ@�i��10.1�j�ɂ��v�Z
      S=0
      do j=1,n
        xj  =a+j*h
        xj_1=a+(j-1)*h
        S=S+f(xj_1)*(xj-xj_1)
      end do
      write(*,901) "S(quadrature by parts 10.1) = ", S

C [�Q�l] �敪���ϖ@�i��10.2�j�ɂ��v�Z
      S=0
      do j=1,n
        xj  =a+j*h
        xj_1=a+(j-1)*h
        S=S+f(xj)*(xj-xj_1)
      end do
      write(*,901) "S(quadrature by parts 10.2) = ", S

  901 format(A30,f10.6)
      end

C**********************************************************************
C �֐�f�̒�`
C**********************************************************************
      real function f(x)
      implicit none
      real :: x

      f=sin(x)
      end function
