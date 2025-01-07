C**********************************************************************
C 3.1�@�x�C���[�@
C**********************************************************************
      implicit none
      integer :: n,j,JJ
      real :: EPS,x,xx,f,ff,fff

C �A���S���Y���菇�P�D
      write(*,*) "Input initial value: (ex.1.0)"
      read(*,*) x
      n=0

      EPS=1.0E-6

      JJ=50 !�ő�50��J��Ԃ��Ă�EPS��菬�����Ȃ�Ȃ�������I��
      do j=1,JJ

C       �A���S���Y���菇�Q�D
        xx=x-f(x)/(ff(x)-f(x)*fff(x)/(2*ff(x))) !��(3.7)�ɕύX

C       �A���S���Y���菇�R�D
        if(abs(xx-x)/abs(x).lt.EPS) then
          exit !�I��
        else
          n=n+1
          x=xx  !�v�Z�����l���R�s�[�i���̌J��Ԃ��̏����j
        end if

C       �i�Q�l�jx�̋����̏o��
         write(*,900) "j=",j,"x=", x

      end do

C     ���ʂ̏o��
  900 format(a3,I2,5x,a3,f10.6)
  901 format(a3,f10.6)
      if(j.ge.JJ) then
        !50��J��Ԃ��Ă�EPS��菬�����Ȃ�Ȃ������ꍇ
        write(*,*) " "
        write(*,*) "It didn't converged."
      else
        write(*,*) " "
        write(*,901) "x=", x
      end if
      end

C**********************************************************************
C �֐��̒�`
C**********************************************************************
      real function f(x)
      real :: x

      f=cos(x)-x**2
      end function
C**********************************************************************
      real function ff(x) !f'(x)�̒�`
      real :: x

      ff=-sin(x)-2*x
      end function
C**********************************************************************
      real function fff(x) !f''(x)�̒�`
      real :: x

      fff=-cos(x)-2
      end function