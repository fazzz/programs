C**********************************************************************
C 11.3　多重積分（一般的な形式）
C**********************************************************************
      implicit none
      integer,parameter :: NM=10
      real :: a,b,c,d,h,k,f,SS
      real :: S(0:NM),x(0:NM)
      integer :: m,n,i,j

      m=10  !x方向にm等分
      n=10  !y方向にn等分
      a=0.0  !x方向の積分区間a～b
      b=1.0
      h=(b-a)/real(m) !x方向刻み幅

      !各XiごとにXiの数値を決める
      do i=0,m
        x(i)=a+i*h
      end do

      !i=0,…,mに対してyに関する積分Sを求める
      do i=0,m
        c=x(i)
        d=2.0*x(i)
        k=(d-c)/real(n)
        S(i)=0.0
        do j=1,n-1
          S(i)=S(i)+f(x(i),c+j*k)
        end do
        S(i)=(f(x(i),c)+2.0*S(i)+f(x(i),d))*k/2.0
      end do

      !Iを計算する
      SS=0.0
      SS=SS+h/2.0*S(0)
      SS=SS+h/2.0*S(m)
	do i=1,m-1
        SS=SS+h*S(i)
      end do

      write(*,901) "I = ", SS

C**********************************************************************
C [参考] 定積分による計算
      write(*,*) " "
      write(*,*) "[for reference]"
      SS=4.0/3.0
      write(*,901) "S(integral from a-b, c-d) = ", SS

  901 format(A30,f15.6)
      end

C**********************************************************************
C 関数fの定義
C**********************************************************************
      real function f(x,y)
      implicit none
      real :: x,y

      f=x+y+1.0
      end function
