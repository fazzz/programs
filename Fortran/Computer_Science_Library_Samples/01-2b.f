C**********************************************************************
C 1.2�@�Q����
C**********************************************************************
      implicit none
      integer :: n,NN
      real :: EPS,x,xx,y,yy,h,gosax,gosay

      NN=700 !�ő��NN��J��Ԃ�
      EPS=1.0E-6 !���Ό덷�������菬�����Ȃ�������������Ƃ���

  900 format(a5,a10,a10)
  901 format(f5.2,f10.6,f10.6) !���ʂ̏o�͂̏���
      write(*,900) "nh", "x", "y"

      x=0.0
      y=1.0
      h=0.01
      do n=1,NN
C       �Q�����̌v�Z
        xx=x+h*y
        yy=y-h*x
        write(*,901) n*h,x,y

C       ��������
        gosax=abs(xx-x)/abs(x)
        gosay=abs(yy-y)/abs(y)
        if(gosax.lt.EPS .and. gosay.lt.EPS) then
          write(*,*) "Converged!"
          exit
        end if

C       ���̎��Ԃ̏���
        x=xx
        y=yy
      end do
      end
