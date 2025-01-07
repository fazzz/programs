C**********************************************************************
C 15.3　ポアソン方程式の差分解法
C**********************************************************************
      implicit none
      integer,parameter :: NM=10
      real :: dx,dy,c
      real :: x(0:NM,0:NM),y(0:NM,0:NM),u(0:NM,0:NM),uu(0:NM,0:NM)
      real :: f(0:NM,0:NM)
      integer :: JJ,J1,KK,K1,j,k,NN,n

C 各定数などの定義
      JJ=10   !x方向の格子数
      J1=JJ-1
      KK=10   !y方向の格子数
      K1=KK-1
      dx=1.0/real(JJ) !領域は1.0×1.0の正方形
      dy=1.0/real(KK)
      c=1.0/(2.0/(dx)**2+2.0/(dy)**2) !係数はまとめておく
      NN=50   !ガウス・ザイデル法の反復回数

C 関数f(x)の定義（配列として保持）
      do k=0,KK
      do j=0,JJ
        f(j,k)=0.0
      end do
      end do

C 初期条件
      do k=0,KK
      do j=0,JJ
        u(j,k)=0.0
       uu(j,k)=0.0
      end do
      end do
      do j=0,JJ
        u(j,0)=1.0
       uu(j,0)=1.0
      end do

C ガウス・ザイデル法の反復
      do n=1,NN

C       反復式の計算
        do k=1,K1
        do j=1,J1
          uu(j,k)=c*( (uu(j-1,k)+u(j+1,k))/(dx)**2
     &               +(uu(j,k-1)+u(j,k+1))/(dy)**2 +f(j,k))
        end do
        end do

C       計算結果のコピー
        do k=1,K1
        do j=1,J1
          u(j,k)=uu(j,k)
        end do
        end do
      end do

C     結果の出力
      do k=0,KK
        do j=0,J1
          write(*,901,advance='no') u(j,k) !改行なし
        end do
        write(*,901) u(JJ,k) !最後だけ改行
      end do
  901 format(f6.3)

      end
