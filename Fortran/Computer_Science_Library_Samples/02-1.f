C**********************************************************************
C 2.1　2分法
C**********************************************************************
      implicit none
      integer :: j,JJ
      real :: a,b,c,x,s,EPS,f

C アルゴリズム手順１．
      a=0.0
      b=1.0
      EPS=1.0E-6

      JJ=50 !最大50回繰り返してもEPSより小さくならなかったら終了
      do j=1,JJ

C       アルゴリズム手順２．
        if(abs(b-a).lt.EPS) then
          x=a  !aまたはbが根
          exit !終了

        else
C         アルゴリズム手順３．
          c=(a+b)/2.0
          s=f(a)*f(c)

C         アルゴリズム手順４．
          if(s.lt.0) then
            b=c
          else if(s.gt.0) then
            a=c
          else
            x=c
            exit
          end if

C       （参考）a,b,cの挙動の出力
         write(*,900) "j=",j,"a=",a,"b=",b,"c=",c
        end if
      end do

C     結果の出力
  900 format(a3,I2,5x,3(a3,f10.6))
  901 format(a3,f10.6)
      if(j.eq.JJ) then
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