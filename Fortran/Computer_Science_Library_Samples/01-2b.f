C**********************************************************************
C 1.2　漸化式
C**********************************************************************
      implicit none
      integer :: n,NN
      real :: EPS,x,xx,y,yy,h,gosax,gosay

      NN=700 !最大でNN回繰り返す
      EPS=1.0E-6 !相対誤差がこれより小さくなったら収束したとする

  900 format(a5,a10,a10)
  901 format(f5.2,f10.6,f10.6) !結果の出力の準備
      write(*,900) "nh", "x", "y"

      x=0.0
      y=1.0
      h=0.01
      do n=1,NN
C       漸化式の計算
        xx=x+h*y
        yy=y-h*x
        write(*,901) n*h,x,y

C       収束判定
        gosax=abs(xx-x)/abs(x)
        gosay=abs(yy-y)/abs(y)
        if(gosax.lt.EPS .and. gosay.lt.EPS) then
          write(*,*) "Converged!"
          exit
        end if

C       次の時間の準備
        x=xx
        y=yy
      end do
      end
