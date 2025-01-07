C**********************************************************************
C 2.4�@�����@
C**********************************************************************
      implicit none
      integer :: n,j,JJ
      real :: EPS,x,xx,xxx,f

C �����l�̎w��
      x =0.0 !x(-1)
      xx=0.5 !x(0)
      n=0

      EPS=1.0E-6

      JJ=50 !�ő�50��J��Ԃ��Ă�EPS��菬�����Ȃ�Ȃ�������I��
      do j=1,JJ

C       �j���[�g���@�̎�(2.4)����(2.6)�ɒu��������
C       x(n+1)��xxx�Cx(n)��xx�Cx(n-1)��x�Ƃ��Ă���
        xxx=xx-(xx-x)/(f(xx)-f(x))*f(xx)

        if(abs(xxx-xx)/abs(xx).lt.EPS) then
          exit !�I��
        else
          n=n+1
          x=xx  !�v�Z�����l���R�s�[�i���̌J��Ԃ��̏����j
          xx=xxx
        end if

C       �i�Q�l�jx�̋����̏o��
         write(*,900) "j=",j,"x=", xx

      end do

C     ���ʂ̕\��
  900 format(a3,I2,5x,a3,f10.6)
  901 format(a3,f10.6)
      if(j.ge.JJ) then
        !50��J��Ԃ��Ă�EPS��菬�����Ȃ�Ȃ������ꍇ
        write(*,*) " "
        write(*,*) "It didn't converged."
      else
        write(*,*) " "
        write(*,901) "x=", xx
      end if
      end

C**********************************************************************
C �֐��̒�`
C**********************************************************************
      real function f(x)
      real :: x

      f=cos(x)-x**2
      end function
