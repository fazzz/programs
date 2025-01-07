C**********************************************************************
C 13.4　高階微分方程式
C**********************************************************************
      implicit none
      integer,parameter :: NM=100
      real :: s1,s2,s3,s4,r1,r2,r3,r4,h,f,g
      real :: x(0:NM),y(0:NM),z(0:NM)
      integer :: n,j

      n=100
      h=0.1
      x(0)=0.0
      y(0)=0.0
      z(0)=0.0

C オイラー法の場合
      do j=0,n-1
        x(j+1)=x(j)+h
        y(j+1)=y(j)+h*f(x(j),y(j),z(j))
        z(j+1)=z(j)+h*g(x(j),y(j),z(j))
      end do
C 結果の出力
  900 format(A4, 3A15)
  901 format(I4, 3f15.6)
      write(*,*) "Eular:"
      write(*,900) "j", "x(j)", "y(j)", "z(j)"
      do j=0,n
        write(*,901) j,x(j),y(j),z(j)
      end do
      write(*,*) " "

C ルンゲ・クッタ法の場合
      do j=0,n-1
        s1=f(x(j),y(j),z(j))
        r1=g(x(j),y(j),z(j))
        s2=f(x(j)+h/2.0,y(j)+h*s1/2.0,z(j)+h*r1/2.0)
        r2=g(x(j)+h/2.0,y(j)+h*s1/2.0,z(j)+h*r1/2.0)
        s3=f(x(j)+h/2.0,y(j)+h*s2/2.0,z(j)+h*r2/2.0)
        r3=g(x(j)+h/2.0,y(j)+h*s2/2.0,z(j)+h*r2/2.0)
        s4=f(x(j)+h,y(j)+h*s3,z(j)+h*r3)
        r4=g(x(j)+h,y(j)+h*s3,z(j)+h*r3)
        x(j+1)=x(j)+h
        y(j+1)=y(j)+h/6.0*(s1+2.0*s2+2.0*s3+s4)
        z(j+1)=z(j)+h/6.0*(r1+2.0*r2+2.0*r3+r4)
      end do
C 結果の出力
      write(*,*) "Runge-Kutta"
      write(*,900) "j", "x(j)", "y(j)", "z(j)"
      do j=0,n
        write(*,901) j,x(j),y(j),z(j)
      end do
      end

C**********************************************************************
C 関数の定義
C**********************************************************************
      real function f(x,y,z)
      real :: x,y,z

      f=z
      end function
C**********************************************************************
      real function g(x,y,z)
      real :: x,y,z

      g=0.25*(1.0-y**2)*z-y+sin(x)
      end function