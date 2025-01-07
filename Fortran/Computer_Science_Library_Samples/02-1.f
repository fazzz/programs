C**********************************************************************
C 2.1�@2���@
C**********************************************************************
      implicit none
      integer :: j,JJ
      real :: a,b,c,x,s,EPS,f

C �A���S���Y���菇�P�D
      a=0.0
      b=1.0
      EPS=1.0E-6

      JJ=50 !�ő�50��J��Ԃ��Ă�EPS��菬�����Ȃ�Ȃ�������I��
      do j=1,JJ

C       �A���S���Y���菇�Q�D
        if(abs(b-a).lt.EPS) then
          x=a  !a�܂���b����
          exit !�I��

        else
C         �A���S���Y���菇�R�D
          c=(a+b)/2.0
          s=f(a)*f(c)

C         �A���S���Y���菇�S�D
          if(s.lt.0) then
            b=c
          else if(s.gt.0) then
            a=c
          else
            x=c
            exit
          end if

C       �i�Q�l�ja,b,c�̋����̏o��
         write(*,900) "j=",j,"a=",a,"b=",b,"c=",c
        end if
      end do

C     ���ʂ̏o��
  900 format(a3,I2,5x,3(a3,f10.6))
  901 format(a3,f10.6)
      if(j.eq.JJ) then
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