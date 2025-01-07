C**********************************************************************
C 10.1　区分求積法と台形公式
C**********************************************************************
      implicit none
      integer :: n,j
      real :: a,b,f,h,S,xj,xj_1

C アルゴリズム手順１．
      n=10
      a=0.0
      b=2.0*atan(1.0) !π/2

C アルゴリズム手順２．
      h=(b-a)/real(n)
      S=0.0

C アルゴリズム手順３．
      do j=1,n-1
        S=S+f(a+j*h)
      end do

C アルゴリズム手順４．
      S=(f(a)+2.0*S+f(b))*h/2.0

C 結果の出力（台形公式の結果）
      write(*,901) "S(trapezoid formula) =        ", S

C**********************************************************************
C [参考] 定積分による計算
      write(*,*) " "
      write(*,*) "[for reference]"
      S=-cos(b)+cos(a)
      write(*,901) "S(integral from a to b) =     ", S

C [参考] 区分求積法（式10.1）による計算
      S=0
      do j=1,n
        xj  =a+j*h
        xj_1=a+(j-1)*h
        S=S+f(xj_1)*(xj-xj_1)
      end do
      write(*,901) "S(quadrature by parts 10.1) = ", S

C [参考] 区分求積法（式10.2）による計算
      S=0
      do j=1,n
        xj  =a+j*h
        xj_1=a+(j-1)*h
        S=S+f(xj)*(xj-xj_1)
      end do
      write(*,901) "S(quadrature by parts 10.2) = ", S

  901 format(A30,f10.6)
      end

C**********************************************************************
C 関数fの定義
C**********************************************************************
      real function f(x)
      implicit none
      real :: x

      f=sin(x)
      end function
