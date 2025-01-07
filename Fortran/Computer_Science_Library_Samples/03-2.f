C**********************************************************************
C 3.2　連立非線形方程式
C**********************************************************************
      implicit none
      integer :: n,j,JJ
      real :: EPS,x,xx,y,yy,dx,dy,aa,f,g,fx,fy,gx,gy

C アルゴリズム手順１．
      write(*,*) "Input initial value for x: (ex.1.0)"
      read(*,*) x
      write(*,*) "Input initial value for y: (ex.0.0)"
      read(*,*) y
      n=0

      EPS=1.0E-6

      JJ=50 !最大50回繰り返してもEPSより小さくならなかったら終了
      do j=1,JJ

C       アルゴリズム手順２．
        f=x**2+y**2-1.0
        g=y-sin(x)
        fx=2.0*x
        fy=2.0*y
        gx=-cos(x)
        gy=1.0

        aa=fy*gx-fx*gy
        dx=(fy*(-g)-gy*(-f))/aa
        dy=(gx*(-f)-fx*(-g))/aa

C       アルゴリズム手順３．
        xx=x+dx
        yy=y+dy

        if(abs(dx).lt.EPS .and. abs(dy).lt.EPS) then
          exit !終了
        else
          n=n+1
          x=xx  !計算した値をコピー（次の繰り返しの準備）
          y=yy
        end if

C       （参考）x,yの挙動の出力
         write(*,900) "j=",j,"x=", x,"y=",y

      end do

C     結果の出力
  900 format(a3,I2,5x,2(a3,f10.6))
  901 format(2(a3,f10.6))
      if(j.ge.JJ) then
        !50回繰り返してもEPSより小さくならなかった場合
        write(*,*) " "
        write(*,*) "It didn't converged."
      else
        write(*,*) " "
        write(*,901) "x=", x, "y=", y
      end if
      end
