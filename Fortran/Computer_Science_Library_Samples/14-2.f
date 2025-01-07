C**********************************************************************
C 14.2　線の方法で1次元熱伝導方程式を解く
C**********************************************************************
      implicit none
      integer,parameter :: NM=10
      real :: dt,dx,c
      real :: x(0:NM),u(0:NM),f(0:NM),uu(0:NM)
      integer :: N,M,j,k

C アルゴリズム手順１．
      N=10  !区間を10等分
      M=100 !t=dt*M=0.2まで計算
      dt=0.002
      dx=1.0/real(N)
      c=dt/(dx)**2 !計算に使う係数をまとめておく

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

C アルゴリズム手順２．
      do j=0,N
        u(j)=f(j)
      end do      

C 結果の出力の準備
      write(*,*) "    t|  x(0)  x(1) ... x(N)"
  901 format(f6.3)

C アルゴリズム手順３．
      do k=1,M
C       ３－１．境界条件
        u(0)=0.0
        u(N)=0.0
C       ３－２．計算
        do j=1,N-1
          uu(j)=u(j)+c*(u(j-1)-2.0*u(j)+u(j+1))
        end do

C       結果の出力
        write(*,901,advance='no') dt*real(k) !t
        write(*,'(a)',advance='no') "| "
        do j=0,N
          write(*,901,advance='no') u(j) !改行なし
        end do
        write(*,*) " " !1行出力終わったら改行

C       uu(j) すなわち n+1 の時の値を今度は u(j) にコピー．
C       次のステップでまた u(j) の値（現時点ではn+1）を用いて
C       uu(j)の値（現時点ではn+2）の値を計算する．
        do j=0,N
          u(j)=uu(j)
        end do
      end do
      end
