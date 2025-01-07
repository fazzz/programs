C**********************************************************************
C 6.3　特殊な大型連立1次方程式（SOR法）
C**********************************************************************
      implicit none
      integer,parameter :: NI=100,NJ=100
      integer :: i,j,II,JJ,m,mmax
      real :: u(NI,NJ)
      real :: gosa,EPS,us,uu,omg

C 問題の設定
      II=100
      JJ=100
      EPS=1.0E-3
      omg=1.2 !0～2の間の値．1だとガウス・ザイデル法に一致

C 初期値
      do j=1,JJ
        do i=1,II
          u(i,j)=0.0
        end do
      end do

C 周囲の温度
      do j=1,JJ
        u( 1,j)=0.5
        u(II,j)=0.0
      end do
      do i=1,II
        u(i, 1)=1.0
        u(i,JJ)=0.0
      end do

C 反復式
      mmax=20000
      do m=1,mmax !反復のためのdo文
        gosa=0.0

        do j=2,JJ-1
          do i=2,II-1
            us=u(i,j) !ν段階の値を記憶（誤差判定に使う）
            uu=(u(i+1,j)+u(i,j+1)+u(i-1,j)+u(i,j-1))/4.0
            gosa=gosa+abs(uu-us)
            u(i,j)=(1.0-omg)*u(i,j)+omg*uu
          end do
        end do

C       収束したかどうかの判定
        if(gosa.lt.EPS) then
          write(*,*) "Converged!"
          write(*,*) " "
          exit !ループから抜ける
        else
          if(mod(m,1000).eq.0) write(*,901) "m=", m, "gosa=", gosa
        end if
      end do

C     mmax回繰り返しても収束しなかった場合
      if(m.gt.mmax) then
        write(*,*) " "
        write(*,*) "It didn't converged."
      end if

C     結果の出力
      open(17,file='06-3b.txt', status='replace')
      do j=1,JJ
        do i=1,II
          write(17,902) 0.01*real(i), 0.01*real(j), u(i,j)
        end do
      end do
      close(17)

  901 format(a2,i6,5x,a5,f10.6)
  902 format(2f5.2,f10.6)
      end

