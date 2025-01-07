C**********************************************************************
C 3.1　ベイリー法
C**********************************************************************
      implicit none
      integer :: n,j,JJ
      real :: EPS,x,xx,f,ff,fff

C アルゴリズム手順１．
      write(*,*) "Input initial value: (ex.1.0)"
      read(*,*) x
      n=0

      EPS=1.0E-6

      JJ=50 !最大50回繰り返してもEPSより小さくならなかったら終了
      do j=1,JJ

C       アルゴリズム手順２．
        xx=x-f(x)/(ff(x)-f(x)*fff(x)/(2*ff(x))) !式(3.7)に変更

C       アルゴリズム手順３．
        if(abs(xx-x)/abs(x).lt.EPS) then
          exit !終了
        else
          n=n+1
          x=xx  !計算した値をコピー（次の繰り返しの準備）
        end if

C       （参考）xの挙動の出力
         write(*,900) "j=",j,"x=", x

      end do

C     結果の出力
  900 format(a3,I2,5x,a3,f10.6)
  901 format(a3,f10.6)
      if(j.ge.JJ) then
        !50回繰り返してもEPSより小さくならなかった場合
        write(*,*) " "
        write(*,*) "It didn't converged."
      else
        write(*,*) " "
        write(*,901) "x=", x
      end if
      end

C**********************************************************************
C 関数の定義
C**********************************************************************
      real function f(x)
      real :: x

      f=cos(x)-x**2
      end function
C**********************************************************************
      real function ff(x) !f'(x)の定義
      real :: x

      ff=-sin(x)-2*x
      end function
C**********************************************************************
      real function fff(x) !f''(x)の定義
      real :: x

      fff=-cos(x)-2
      end function