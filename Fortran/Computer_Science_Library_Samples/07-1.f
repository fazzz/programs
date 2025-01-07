C**********************************************************************
C 7.1　ベキ乗法
C**********************************************************************
      implicit none
      integer,parameter :: NM=3
      integer :: i,j,k,kmax,n
      real :: a(NM,NM),x(NM),y(NM)
      real :: gosa,EPS,lambda,lambda2,sum1,sum2

      EPS=1.0E-6

C アルゴリズム手順１．
      n=3
      a(1,1)=11.0; a(1,2)= 7.0; a(1,3)=-5.0
      a(2,1)= 0.0; a(2,2)=10.0; a(2,3)=-1.0
      a(3,1)= 2.0; a(3,2)= 8.0; a(3,3)= 3.0
 
      !適当でよいが、ゼロベクトルにすると計算できない
      do i=1,n
        x(i)=1.0
      end do

C アルゴリズム手順２．
      kmax=200
      lambda=10000.0 !λの初期値（ダミー）
      lambda2=10000.0
      do k=1,kmax
        gosa=0.0

C       (1) y=Ax
        do i=1,n
          y(i)=0.0
          do j=1,n
            y(i)=y(i)+a(i,j)*x(j)
          end do
        end do

C       (1) λ=(y,y)/(y,x)
        sum1=0.0
        sum2=0.0
        do i=1,n
          sum1=sum1+y(i)*y(i)
          sum2=sum2+y(i)*x(i)
        end do
        lambda=sum1/sum2

C       (2)
        do i=1,n
          x(i)=y(i)/sqrt(sum1)
        end do

C       アルゴリズム手順３．
        gosa=abs(lambda-lambda2)/abs(lambda2)
        lambda2=lambda
        if(gosa.lt.EPS) then
          write(*,*) "Converged!"
          write(*,*) " "
          exit !ループから抜ける
        end if

C       結果の出力
        write(*,901) k,lambda
      end do

C     kmax回繰り返しても収束しなかった場合
      if(k.gt.kmax) then
        write(*,*) " "
        write(*,*) "It didn't converged."
      end if

  901 format(I3,f10.6)
      end

