C**********************************************************************
C 2.4　割線法
C**********************************************************************
      implicit none
      integer :: n,j,JJ
      real :: EPS,x,xx,xxx,f

C 初期値の指定
      x =0.0 !x(-1)
      xx=0.5 !x(0)
      n=0

      EPS=1.0E-6

      JJ=50 !最大50回繰り返してもEPSより小さくならなかったら終了
      do j=1,JJ

C       ニュートン法の式(2.4)を式(2.6)に置き換える
C       x(n+1)をxxx，x(n)をxx，x(n-1)をxとしている
        xxx=xx-(xx-x)/(f(xx)-f(x))*f(xx)

        if(abs(xxx-xx)/abs(xx).lt.EPS) then
          exit !終了
        else
          n=n+1
          x=xx  !計算した値をコピー（次の繰り返しの準備）
          xx=xxx
        end if

C       （参考）xの挙動の出力
         write(*,900) "j=",j,"x=", xx

      end do

C     結果の表示
  900 format(a3,I2,5x,a3,f10.6)
  901 format(a3,f10.6)
      if(j.ge.JJ) then
        !50回繰り返してもEPSより小さくならなかった場合
        write(*,*) " "
        write(*,*) "It didn't converged."
      else
        write(*,*) " "
        write(*,901) "x=", xx
      end if
      end

C**********************************************************************
C 関数の定義
C**********************************************************************
      real function f(x)
      real :: x

      f=cos(x)-x**2
      end function
