C**********************************************************************
C 11.3　多重積分（台形公式）
C**********************************************************************
      implicit none
      integer :: i,j,m,n
      real :: a,b,c,d,h,k,S,f,PAI

      PAI=4.0*atan(1.0)
      m=40  !x方向にm等分
      n=30  !y方向にn等分
      a=-2.0  !x方向の積分区間a～b
      b=2.0
      c=0.0   !y方向の積分区間c～d
      d=PAI
      h=(b-a)/real(m) !x方向刻み幅
      k=(d-c)/real(n) !y方向刻み幅

      S=0.0

      !積分領域内
      do j=1,n-1
        do i=1,m-1
          S=S+4.0*f(a+i*h,c+j*k)
        end do
      end do

      !4頂点を除く境界
      do i=1,m-1
        S=S+2.0*f(a+i*h,c)
      end do
      do i=1,m-1
        S=S+2.0*f(a+i*h,d)
      end do
      do j=1,n-1
        S=S+2.0*f(a,c+j*k)
      end do
      do j=1,n-1
        S=S+2.0*f(b,c+j*k)
      end do

      !境界の4頂点
      S=S+1.0*f(a,c)
      S=S+1.0*f(a,d)
      S=S+1.0*f(b,c)
      S=S+1.0*f(b,d)

      S=h*k/4.0*S
      write(*,901) "S(trapezoid formula) =        ", S

C**********************************************************************
C [参考] 定積分による計算
      write(*,*) " "
      write(*,*) "[for reference]"
      S=88.0/3.0*PAI
      write(*,901) "S(integral from a-b, c-d) =   ", S

  901 format(A30,f15.6)
      end

C**********************************************************************
C 関数fの定義
C**********************************************************************
      real function f(x,y)
      implicit none
      real :: x,y

      f=(x+1)**2-cos(3*y)+5.0
      end function
