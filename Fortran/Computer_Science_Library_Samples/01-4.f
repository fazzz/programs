C**********************************************************************
C 1.4　桁落ちと情報落ち
C**********************************************************************
      implicit none
      real :: a,b,c,D,x1
      double precision :: aa,bb,cc,DD,xx1

      do while(.true.)
C アルゴリズム手順１．
        write(*,*) "Input a,b,c (axx+bx+c)"
        write(*,*) " ex1. a=-1.0, b=-1.0, c= 6.0"
        write(*,*) " ex2. a= 2.0, b= 4.0, c= 2.0"
        write(*,*) " ex3. a= 2.0, b= 2.0, c= 1.0"
        write(*,*) "[for reference] a=1.0, b=100000.0, c=2.0"
         read(*,*) a,b,c

C アルゴリズム手順２．
        if(a.eq.0.0) then
          write(*,*) "Input again (""a"" shouldn't be 0.0)"
          write(*,*) " "
          continue
        else
          exit
        end if
      end do

C アルゴリズム手順３．
      D=b**2.0-4.0*a*c

C     アルゴリズム手順４．
      if(D.lt.0.0) then
        write(*,*) "two imaginary solutions"
        write(*,*) " real part     :", -b/(2.0*a)
        write(*,*) " imaginary part:", sqrt(abs(D))/(2.0*a)

C     アルゴリズム手順５．
      else if(D.eq.0.0) then
        write(*,*) "multiple root"
        write(*,*) "x=", -b/(2.0*a)

      else
C       アルゴリズム手順６．
        if(b.ge.0.0) then
          write(*,*) "two real solutions"
          x1=(-b-sqrt(D))/(2.0*a)
          write(*,*) "x1=", x1
          write(*,*) "x2=", c/(a*x1)

C       アルゴリズム手順７．
        else
          write(*,*) "two real solutions"
          x1=(-b+sqrt(D))/(2.0*a)
          write(*,*) "x1=", x1
          write(*,*) "x2=", c/(a*x1)
        end if
      end if

C**********************************************************************
C （参考）bの絶対値が大きい場合、桁落ち・情報落ちが起こることがある
C**********************************************************************
      if(D.gt.0.0 .and. abs(b).gt.abs(a*c*100.0)) then

C 2次方程式の根の公式を利用した計算
        write(*,*) " "
        write(*,*) "[for reference]"
        write(*,*) " fear of error"
        if(b.ge.0.0) then
          write(*,*) "  x1=",(-b-sqrt(D))/(2.0*a)
          write(*,*) "  x2=",(-b+sqrt(D))/(2.0*a)
        else
          write(*,*) "  x1=",(-b+sqrt(D))/(2.0*a)
          write(*,*) "  x2=",(-b-sqrt(D))/(2.0*a)
        end if

C 倍精度実数での計算
        aa=dble(a)
        bb=dble(b)
        cc=dble(c)
        DD=bb**2-4*aa*cc

        write(*,*) " solution by double precision"
        if(b.ge.0.0) then
          write(*,*) "  x1=",(-bb-sqrt(DD))/(2.0*aa)
          write(*,*) "  x2=",(-bb+sqrt(DD))/(2.0*aa)
        else
          write(*,*) "  x1=",(-bb+sqrt(DD))/(2.0*aa)
          write(*,*) "  x2=",(-bb-sqrt(DD))/(2.0*aa)
        end if
      end if
      end
