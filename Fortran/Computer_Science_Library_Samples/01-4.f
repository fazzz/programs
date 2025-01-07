C**********************************************************************
C 1.4�@�������Ə�񗎂�
C**********************************************************************
      implicit none
      real :: a,b,c,D,x1
      double precision :: aa,bb,cc,DD,xx1

      do while(.true.)
C �A���S���Y���菇�P�D
        write(*,*) "Input a,b,c (axx+bx+c)"
        write(*,*) " ex1. a=-1.0, b=-1.0, c= 6.0"
        write(*,*) " ex2. a= 2.0, b= 4.0, c= 2.0"
        write(*,*) " ex3. a= 2.0, b= 2.0, c= 1.0"
        write(*,*) "[for reference] a=1.0, b=100000.0, c=2.0"
         read(*,*) a,b,c

C �A���S���Y���菇�Q�D
        if(a.eq.0.0) then
          write(*,*) "Input again (""a"" shouldn't be 0.0)"
          write(*,*) " "
          continue
        else
          exit
        end if
      end do

C �A���S���Y���菇�R�D
      D=b**2.0-4.0*a*c

C     �A���S���Y���菇�S�D
      if(D.lt.0.0) then
        write(*,*) "two imaginary solutions"
        write(*,*) " real part     :", -b/(2.0*a)
        write(*,*) " imaginary part:", sqrt(abs(D))/(2.0*a)

C     �A���S���Y���菇�T�D
      else if(D.eq.0.0) then
        write(*,*) "multiple root"
        write(*,*) "x=", -b/(2.0*a)

      else
C       �A���S���Y���菇�U�D
        if(b.ge.0.0) then
          write(*,*) "two real solutions"
          x1=(-b-sqrt(D))/(2.0*a)
          write(*,*) "x1=", x1
          write(*,*) "x2=", c/(a*x1)

C       �A���S���Y���菇�V�D
        else
          write(*,*) "two real solutions"
          x1=(-b+sqrt(D))/(2.0*a)
          write(*,*) "x1=", x1
          write(*,*) "x2=", c/(a*x1)
        end if
      end if

C**********************************************************************
C �i�Q�l�jb�̐�Βl���傫���ꍇ�A�������E��񗎂����N���邱�Ƃ�����
C**********************************************************************
      if(D.gt.0.0 .and. abs(b).gt.abs(a*c*100.0)) then

C 2���������̍��̌����𗘗p�����v�Z
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

C �{���x�����ł̌v�Z
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
