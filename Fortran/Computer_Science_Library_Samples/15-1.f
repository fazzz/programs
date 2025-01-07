C**********************************************************************
C 15.1　移流方程式の差分解法
C**********************************************************************
      implicit none
      integer,parameter :: NM=20
      real :: dt,dx,c,r
      real :: x(0:NM),u(0:NM),f(0:NM),uu(0:NM)
      integer :: JJ,NN,j,n

C 各定数などの定義
      JJ=20            !x方向の離散点の数
      dx=2.0/real(JJ)  !領域をJJ等分
      dt=0.01          !dtは指定
      NN=30            !t方向の離散点の数（NNステップ目まで計算）
      c=2.0            !移流速度
      r=c*dt/dx

C 離散点の定義
      do j=0,JJ
        x(j)=real(j)*dx
      end do

C 関数f(x)の定義（配列として保持）
      do j=0,JJ
        f(j)=sin(3.0*x(j))
        if(f(j).lt.0.0) f(j)=0.0
      end do

C 初期条件（式(15.3)）
      do j=0,JJ
        u(j)=f(j)
       uu(j)=f(j) 
      end do      
      call out(dt,n,JJ,u) !初期条件の出力

C 時間に対して1からNNまで1行ずつ計算
      do n=1,NN
C       式(15.4)の計算
C       左辺（t=n+1のときの値）を配列uuに格納する
        do j=1,JJ
          uu(j)=(1-r)*u(j)+r*u(j-1)
        end do

C       計算結果uuの出力
        call out(dt,n,JJ,uu)

C       uu(j) すなわち n+1 の時の値を u(j) にコピー．
C       （次の時間の準備）
        do j=0,JJ
          u(j)=uu(j)
        end do
      end do !これで時間の計算が一行終了
      end

C**********************************************************************
C 結果を出力するためのサブルーチン
C**********************************************************************
      subroutine out(dt,n,JJ,u)
      real :: dt
      real :: u(JJ)
      integer :: n,j,JJ

      write(*,901,advance='no') dt*real(n) !t
      write(*,'(a)',advance='no') "| "
      do j=0,JJ-1
        write(*,902,advance='no') u(j) !改行なし
      end do
      write(*,902) u(JJ) !最後だけ改行

  901 format(f6.2)
  902 format(f6.3)
      end subroutine
