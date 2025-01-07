C**********************************************************************
C 6.2　SOR法
C**********************************************************************
      implicit none
      integer,parameter :: NM=4
      integer :: j,n,k,m,mmax
      real :: a(NM,NM),b(NM),x(NM)
      real :: gosa,EPS,xs,omg

C 問題の設定
      n=4
      a(1,1)= 9.0; a(1,2)= 2.0; a(1,3)= 1.0; a(1,4)= 1.0;  b(1)=20.0
      a(2,1)= 2.0; a(2,2)= 8.0; a(2,3)=-2.0; a(2,4)= 1.0;  b(2)=16.0
      a(3,1)=-1.0; a(3,2)=-2.0; a(3,3)= 7.0; a(3,4)=-2.0;  b(3)= 8.0
      a(4,1)= 1.0; a(4,2)=-1.0; a(4,3)=-2.0; a(4,4)= 6.0;  b(4)=17.0
 
      EPS=1.0E-6
      omg=1.2 !0～2の間の値．1だとガウス・ザイデル法に一致

C 初期値
      do j=1,n
        x(j)=0.0
      end do

C 反復式
      mmax=50
      do m=1,mmax !反復のためのdo文
        gosa=0.0

        do j=1,n  !各Xjに対するdo文
          xs=b(j) !仮の値
          do k=1,j-1
            xs=xs-a(j,k)*x(k)
          end do
          do k=j+1,n
            xs=xs-a(j,k)*x(k)
          end do
          xs=xs/a(j,j)

          gosa=gosa+abs(xs-x(j))/abs(x(j))
          x(j)=(1.0-omg)*x(j)+omg*xs

        end do

C       結果の出力
        write(*,901,advance='no') "Nu=",m+1
        do j=1,n
          write(*,902,advance='no') "x(",j,")=",x(j)
        end do
        write(*,*) " "

C       収束したかどうかの判定
        if(gosa.lt.EPS) then
          write(*,*) "Converged!"
          write(*,*) " "
          exit !ループから抜ける
        end if
      end do

C     mmax回繰り返しても収束しなかった場合
      if(m.gt.mmax) then
        write(*,*) " "
        write(*,*) "It didn't converged."
      end if

  901 format(a3,I2)
  902 format(a3,I2,a2,f10.6)
      end

