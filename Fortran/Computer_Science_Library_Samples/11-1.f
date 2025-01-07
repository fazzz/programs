C**********************************************************************
C 11.1�@�����o�[�O�ϕ�
C**********************************************************************
      implicit none
      integer,parameter :: NM=10
      integer :: n,m,j,i
      real :: a,b,f,h,Stmp,SS
      real :: S(0:NM,0:NM)

      n=50
      m=5 !�ŏI�I�ɋ�Ԃ� 2**m * n��������
      a=0.0
      b=5.0


! S(0,0), S(1,0), �c, S(m,0)���܂��v�Z����
      h=(b-a)/2.0**0
      S(0,0)=h/2.0*(f(a)+f(b))
      do j=1,m
        h=(b-a)/2.0**j
        Stmp=0.0
        do i=1,2**(j-1)
          Stmp=Stmp+f(a+(2*i-1)*h)
        end do
        S(j,0)=1.0/2.0*S(j-1,0) + h*Stmp
      end do

! S(m,j)�����Ɍv�Z����
      do j=1,m
        do n=j,m
          S(n,j)=(4.0**j*S(n,j-1)-S(n-1,j-1))/(4.0**j-1.0)
        end do
      end do

! ���ʂ̏o��
      write(*,901) "S(Romberg integration) =       ", S(m,m)

C**********************************************************************
C [�Q�l] �s��ϕ��ɂ��v�Z
      write(*,*) " "
      write(*,*) "[for reference]"
      SS=(sin(b)+1.0/2.0*b**2)-(sin(a)+1.0/2.0*a**2)
      write(*,901) "S(integral from a to b) =     ", SS

C [�Q�l] ��`�����ɂ��v�Z
      h=(b-a)/real(n)
      SS=0.0
      do j=1,n-1
        SS=SS+f(a+j*h)
      end do
      SS=(f(a)+2.0*SS+f(b))*h/2.0
      write(*,901) "S(trapezoid formula) =        ", SS

  901 format(A30,f10.6)
      end

C**********************************************************************
C �֐�f�̒�`
C**********************************************************************
      real function f(x)
      implicit none
      real :: x

      f=cos(x)+x
      end function
