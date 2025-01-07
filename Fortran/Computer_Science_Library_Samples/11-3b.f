C**********************************************************************
C 11.3�@���d�ϕ��i��ʓI�Ȍ`���j
C**********************************************************************
      implicit none
      integer,parameter :: NM=10
      real :: a,b,c,d,h,k,f,SS
      real :: S(0:NM),x(0:NM)
      integer :: m,n,i,j

      m=10  !x������m����
      n=10  !y������n����
      a=0.0  !x�����̐ϕ����a�`b
      b=1.0
      h=(b-a)/real(m) !x�������ݕ�

      !�eXi���Ƃ�Xi�̐��l�����߂�
      do i=0,m
        x(i)=a+i*h
      end do

      !i=0,�c,m�ɑ΂���y�Ɋւ���ϕ�S�����߂�
      do i=0,m
        c=x(i)
        d=2.0*x(i)
        k=(d-c)/real(n)
        S(i)=0.0
        do j=1,n-1
          S(i)=S(i)+f(x(i),c+j*k)
        end do
        S(i)=(f(x(i),c)+2.0*S(i)+f(x(i),d))*k/2.0
      end do

      !I���v�Z����
      SS=0.0
      SS=SS+h/2.0*S(0)
      SS=SS+h/2.0*S(m)
	do i=1,m-1
        SS=SS+h*S(i)
      end do

      write(*,901) "I = ", SS

C**********************************************************************
C [�Q�l] ��ϕ��ɂ��v�Z
      write(*,*) " "
      write(*,*) "[for reference]"
      SS=4.0/3.0
      write(*,901) "S(integral from a-b, c-d) = ", SS

  901 format(A30,f15.6)
      end

C**********************************************************************
C �֐�f�̒�`
C**********************************************************************
      real function f(x,y)
      implicit none
      real :: x,y

      f=x+y+1.0
      end function
