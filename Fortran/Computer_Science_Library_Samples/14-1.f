C**********************************************************************
C 14.1　2階微分方程式の境界値問題
C**********************************************************************
      implicit none
      integer,parameter :: NM=20
      real :: p,q,r,AA,BB,h,exact
      real :: x(0:NM),y(0:NM)
      real :: a(0:NM),b(0:NM),c(0:NM),d(0:NM)
      integer :: JJ,J1,j

C アルゴリズム手順１．２．その他宣言
      AA=0.0 !アルゴリズム中の「A」
      BB=0.0 !アルゴリズム中の「B」
      JJ=20  !アルゴリズム中の「J」
      x( 0)=0.0
      x(JJ)=1.0
      y( 0)=AA
      y(JJ)=BB
      h=(x(JJ)-x(0))/real(JJ)
      J1=JJ-1

C アルゴリズム手順３．
      do j=1,J1
        x(j)=x(0)+real(j)*h
        a(j)=1.0-h*p(x(j))/2.0
        b(j)=-(2.0-h**2*q(x(j)))
        c(j)=1.0+h*p(x(j))/2.0
        d(j)=h**2*r(x(j))
      end do

C アルゴリズム手順４．
      d(1)=d(1)-a(1)*AA
      d(J1)=d(J1)-c(J1)*BB

C アルゴリズム手順５．
      call THOMAS(1,J1,a,b,c,d,y)

C 結果の出力
      write(*,900) "x(j)", "calculated", "exact"
  900 format(A6, 2A15)
      do j=0,JJ
        write(*,901) x(j),y(j),exact(x(j))
      end do
  901 format(f6.2, 2f15.6)
      end

C**********************************************************************
C トーマス法を計算するためのサブルーチン
C**********************************************************************
      subroutine THOMAS(IL,IU,a,b,c,d,y)
      implicit none
      integer,parameter :: NM=20
      integer :: IL,IU,j
      real :: a(0:NM),b(0:NM),c(0:NM),d(0:NM)
      real :: g(0:NM),s(0:NM)
      real :: x(0:NM),y(0:NM)

      g(IL)=b(IL)
      s(IL)=d(IL)
      do j=IL+1,IU
        g(j)=b(j)-a(j)*c(j-1)/g(j-1)
        s(j)=d(j)-a(j)*s(j-1)/g(j-1)
      end do
      y(IU)=s(IU)/g(IU)
      do j=IU-1,IL,-1
        y(j)=(s(j)-c(j)*y(j+1))/g(j)
      end do
      end subroutine

C**********************************************************************
C 関数の定義
C**********************************************************************
      real function p(x)
      implicit none
      real :: x
      p=2.0
      end function
C**********************************************************************
      real function q(x)
      implicit none
      real :: x
      q=2.0
      end function
C**********************************************************************
      real function r(x)
      implicit none
      real :: x
      r=4.0*x
      end function
C**********************************************************************
C 厳密解を計算するためのサブルーチン
C**********************************************************************
      real function exact(x)
      real :: x

      exact=2.0*exp(-x)*(cos(x)-sin(x)*cos(1.0)/sin(1.0))+2.0*x-2.0
      end function
