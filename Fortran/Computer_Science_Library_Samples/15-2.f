C**********************************************************************
C 15.2　拡散方程式の差分解法
C**********************************************************************
      implicit none
      integer,parameter :: NM=10
      real :: dt,dx,r
      real :: x(0:NM),u(0:NM),f(0:NM),uu(0:NM)
      real :: a(0:NM),b(0:NM),c(0:NM),d(0:NM)
      integer :: N,N1,M,j,k

C 各定数などの定義
      N=10
      N1=N-1
      M=100
      dt=0.002
      dx=1.0/real(N)
      r=dt/(dx)**2

C 離散点の定義
      do j=0,N
        x(j)=real(j)*dx
      end do

C 関数f(x)の定義（配列として保持）
      do j=0,N
        if(x(j).le.0.5) then
          f(j)=2.0*x(j)
        else
          f(j)=2.0*(1.0-x(j))
        end if
      end do

C 初期条件
      do j=0,N
        u(j)=f(j)
       uu(j)=f(j)
      end do
      call out(dt,0,N,u)

C 3項方程式の係数を定義（a,b,cは定数。dは時間ごとに更新）
      do j=1,N1
        a(j)=-r
        b(j)=1.0+2.0*r
        c(j)=-r
      end do

C 時間に対する計算
      do k=1,M

C       3項方程式の係数を更新
        do j=1,N1
          d(j)=u(j)
        end do
        d( 1)=d( 1)-a( 1)*u(0)
        d(N1)=d(N1)-c(N1)*u(N)

C       トーマス法
        call THOMAS(1,N1,a,b,c,d,uu) !結果を配列uuに格納

C       結果の出力
        call out(dt,k,N,uu)

C       uuをuにコピー（次の時間の準備）
        do j=1,N1
          u(j)=uu(j)
        end do
      end do
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

C**********************************************************************
C 結果を出力するためのサブルーチン
C**********************************************************************
      subroutine out(dt,n,JJ,u)
      real :: dt
      real :: u(0:JJ)
      integer :: n,j,JJ

      write(*,901,advance='no') dt*real(n) !t
      write(*,'(a)',advance='no') "| "
      do j=0,JJ-1
        write(*,901,advance='no') u(j) !改行なし
      end do
      write(*,901) u(JJ) !最後だけ改行

  901 format(f6.3)
      end subroutine
