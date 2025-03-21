C**********************************************************************
C 11.3@½dÏªiä`ö®j
C**********************************************************************
      implicit none
      integer :: i,j,m,n
      real :: a,b,c,d,h,k,S,f,PAI

      PAI=4.0*atan(1.0)
      m=40  !xûüÉmª
      n=30  !yûüÉnª
      a=-2.0  !xûüÌÏªæÔa`b
      b=2.0
      c=0.0   !yûüÌÏªæÔc`d
      d=PAI
      h=(b-a)/real(m) !xûüÝ
      k=(d-c)/real(n) !yûüÝ

      S=0.0

      !ÏªÌæà
      do j=1,n-1
        do i=1,m-1
          S=S+4.0*f(a+i*h,c+j*k)
        end do
      end do

      !4¸_ð­«E
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

      !«EÌ4¸_
      S=S+1.0*f(a,c)
      S=S+1.0*f(a,d)
      S=S+1.0*f(b,c)
      S=S+1.0*f(b,d)

      S=h*k/4.0*S
      write(*,901) "S(trapezoid formula) =        ", S

C**********************************************************************
C [Ql] èÏªÉæévZ
      write(*,*) " "
      write(*,*) "[for reference]"
      S=88.0/3.0*PAI
      write(*,901) "S(integral from a-b, c-d) =   ", S

  901 format(A30,f15.6)
      end

C**********************************************************************
C ÖfÌè`
C**********************************************************************
      real function f(x,y)
      implicit none
      real :: x,y

      f=(x+1)**2-cos(3*y)+5.0
      end function
