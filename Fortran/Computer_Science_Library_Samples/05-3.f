C**********************************************************************
C 5.3　トーマス法
C**********************************************************************
      implicit none
      integer,parameter :: NM=5
      integer :: i,n
      real :: a(NM),b(NM),c(NM),d(NM)
      real :: g(NM),s(NM),x(NM)

C 問題の設定
      n=5
      do i=1,n
        a(i)=1.0
        b(i)=2.0
        c(i)=1.0
        d(i)=real(i)
      end do

C     アルゴリズム手順１．
      g(1)=b(1)
      s(1)=d(1)

C     アルゴリズム手順２．
      do i=2,n
        g(i)=b(i)-a(i)*c(i-1)/g(i-1)
        s(i)=d(i)-a(i)*s(i-1)/g(i-1)
      end do

C     アルゴリズム手順３．
      x(n)=s(n)/g(n)

C     アルゴリズム手順４．
      do i=n-1,1,-1
        x(i)=(s(i)-c(i)*x(i+1))/g(i)
      end do

C 結果の出力
      do i=1,n
        write(*,901) "x(",i,")=",x(i)
      end do
  901 format(A2,I2,A2,f15.6)
      end

