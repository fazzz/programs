C**********************************************************************
C 9.1　スプライン補間
C**********************************************************************
      implicit none
      integer,parameter :: NM=10
      integer,parameter :: NA=20
      integer :: n,k,i,j
      real :: x(0:NM),f(0:NM),xx(0:NA)
      real :: a(0:NM),b(0:NM),c(0:NM),d(0:NM)
      real :: S2(0:NM)
      real :: SS,Hk

      write(*,902) "xx", "S(xx)"
  902 format(a3,a10)

C アルゴリズム手順１．
      !値が分かっている離散点たち
      n=8  !区間数
      do i=0,n
        x(i)=1.0+0.25*real(i)
        f(i)=1.0/(1.0+10.0*x(i)**8)
      end do

      !1から3まで0.1刻みで21点分の「xx」に対して「S(xx)」を計算する
      do i=0,NA
        xx(i)=1.0+0.1*real(i)
      end do

C アルゴリズム手順２．
      !Sの2階微分に対する連立1次方程式を立てる
      !3項方程式になるのでトーマス法を利用
      do k=1,n-1
        a(k)=x(k)-x(k-1)
        b(k)=2.0*(x(k+1)-x(k-1))
        c(k)=x(k+1)-x(k)
        d(k)=6.0*( (f(k+1)-f(k))/(x(k+1)-x(k))
     &            -(f(k)-f(k-1))/(x(k)-x(k-1)))
      end do
      a(1)=0.0
      c(n-1)=0.0

C アルゴリズム手順３．
      S2(0)=0.0
      S2(n)=0.0

C アルゴリズム手順４．
      call THOMAS(1,n-1,a,b,c,d,S2)

C 各xx(i)に対して、どの小区間にいるのか
C すなわち、どの補間式S(k)を使うべきか探す
      do i=0,NA
        if(xx(i).lt.x(0)) then
          k=0
        else if(x(n).le.xx(i)) then
          k=n-1
        else
          do j=0,n-1
          if(x(j).le.xx(i) .and. xx(i).lt.x(j+1)) then
            k=j
            exit
          end if
          end do
        end if

C アルゴリズム手順５．
        Hk=x(k+1)-x(k)
        SS= (S2(k+1)*(xx(i)-x(k))**3 - S2(k)*(xx(i)-x(k+1))**3)/(6.0*Hk)
     &     +(f(k+1)/Hk - Hk*S2(k+1)/6.0)*(xx(i)-x(k))
     &     -(f(k)/Hk - Hk*S2(k)/6.0)*(xx(i)-x(k+1))
        write(*,901) xx(i),SS
      end do
  901 format(f3.1,f10.6)

      end

C**********************************************************************
C トーマス法を計算するためのサブルーチン
C**********************************************************************
      subroutine THOMAS(IL,IU,a,b,c,d,y)
      implicit none
      integer,parameter :: NM=10
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
