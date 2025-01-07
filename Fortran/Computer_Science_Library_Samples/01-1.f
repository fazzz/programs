C**********************************************************************
C 1.1�@�A���S���Y��
C**********************************************************************
      implicit none
      integer,parameter :: NM=10
      real :: a(0:NM),y(0:NM),xx(0:NM)
      real :: x,P
      integer :: j,n,na,nt

C �W��a�̎w��
      a( 0)= 0.3; a( 1)= 0.2; a( 2)=-1.0; a( 3)=-0.5; a( 4)= 1.0
      a( 5)=-8.0; a( 6)=-0.1; a( 7)= 9.0; a( 8)= 1.1; a( 9)=12.3
      a(10)= 2.4
      n=10

C �A���S���Y���菇�P�D
      write(*,*) "Input x:"
       read(*,*) x

C �A���S���Y���菇�Q�D
      y(0)=a(0)

C �A���S���Y���菇�R�D
      nt=0
      na=0
      do j=1,n
        y(j)=y(j-1)*x+a(j)
        nt=nt+1 !�i�Q�l�j��Z�̉񐔂��L��
        na=na+1 !�i�Q�l�j���Z�̉񐔂��L��
      end do

C ���ʂ̏o��
      write(*,*) "P(x)=", y(n)
      write(*,*) "Multiplication: ", nt
      write(*,*) "Addition:       ", na

C**********************************************************************
C�i�Q�l�j��(1.1)�����̂܂܌v�Z�����ꍇ
C**********************************************************************
      nt=0
      na=0
      xx(1)=x
      do j=2,n
        xx(j)=xx(j-1)*x
        nt=nt+1
      end do

      P=a(n)
      do j=n-1,0,-1
        P=P+a(j)*xx(n-j)
        nt=nt+1
        na=na+1
      end do

C ���ʂ̏o��
      write(*,*) " "
      write(*,*) "[for reference]"
      write(*,*) "P(x)=", P
      write(*,*) "Multiplication: ", nt
      write(*,*) "Addition:       ", na

      end
