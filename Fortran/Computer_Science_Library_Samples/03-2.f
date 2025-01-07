C**********************************************************************
C 3.2�@�A������`������
C**********************************************************************
      implicit none
      integer :: n,j,JJ
      real :: EPS,x,xx,y,yy,dx,dy,aa,f,g,fx,fy,gx,gy

C �A���S���Y���菇�P�D
      write(*,*) "Input initial value for x: (ex.1.0)"
      read(*,*) x
      write(*,*) "Input initial value for y: (ex.0.0)"
      read(*,*) y
      n=0

      EPS=1.0E-6

      JJ=50 !�ő�50��J��Ԃ��Ă�EPS��菬�����Ȃ�Ȃ�������I��
      do j=1,JJ

C       �A���S���Y���菇�Q�D
        f=x**2+y**2-1.0
        g=y-sin(x)
        fx=2.0*x
        fy=2.0*y
        gx=-cos(x)
        gy=1.0

        aa=fy*gx-fx*gy
        dx=(fy*(-g)-gy*(-f))/aa
        dy=(gx*(-f)-fx*(-g))/aa

C       �A���S���Y���菇�R�D
        xx=x+dx
        yy=y+dy

        if(abs(dx).lt.EPS .and. abs(dy).lt.EPS) then
          exit !�I��
        else
          n=n+1
          x=xx  !�v�Z�����l���R�s�[�i���̌J��Ԃ��̏����j
          y=yy
        end if

C       �i�Q�l�jx,y�̋����̏o��
         write(*,900) "j=",j,"x=", x,"y=",y

      end do

C     ���ʂ̏o��
  900 format(a3,I2,5x,2(a3,f10.6))
  901 format(2(a3,f10.6))
      if(j.ge.JJ) then
        !50��J��Ԃ��Ă�EPS��菬�����Ȃ�Ȃ������ꍇ
        write(*,*) " "
        write(*,*) "It didn't converged."
      else
        write(*,*) " "
        write(*,901) "x=", x, "y=", y
      end if
      end
