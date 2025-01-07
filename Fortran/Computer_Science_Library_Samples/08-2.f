C**********************************************************************
C 8.2�@�G���~�[�g���
C**********************************************************************
      implicit none
      integer,parameter :: NM=2
      integer :: n,j,jj
      real :: x(0:NM),f(0:NM),ff(0:NM)
      real :: l(0:NM),d(0:NM),h(0:NM),g(0:NM)
      real :: L_bunshi1,L_bunshi2,L_bunbo1,L_bunbo2,HH,a

C �A���S���Y���菇�P�D
      n=NM  !n+1�_��ʂ�2n+1������p����

      !�l���������Ă��闣�U�_����
      x(0)=0.0
      f(0)=0.000000
      ff(0)=1.000000
      x(1)=0.1
      f(1)=0.099833
      ff(1)=0.995004
      x(2)=0.2
      f(2)=0.198669
      ff(2)=0.980067

      a=0.15

C     �A���S���Y���菇�Q�D
      do j=0,n
        L_bunshi1=1.0
        L_bunshi2=1.0
        L_bunbo1 =1.0
        L_bunbo2 =1.0
        do jj=0,j-1
          L_bunshi1=L_bunshi1*(a   -x(jj))   !(a-X0)(a-X1)�c(a-Xj-1)
          L_bunbo1 =L_bunbo1 *(x(j)-x(jj))   !(Xj-X0)(Xj-X1)�c(Xj-Xj-1)
        end do
        do jj=j+1,n
          L_bunshi2=L_bunshi2*(a   -x(jj))   !(a-Xj+1)�c(a-Xn)
          L_bunbo2 =L_bunbo2 *(x(j)-x(jj))   !(Xj-Xj+1)�c(Xj-Xn)
        end do
        L(j)=(L_bunshi1*L_bunshi2)/(L_bunbo1*L_bunbo2)
      end do

      do j=0,n
        d(j)=0.0
        do jj=0,j-1
          d(j)=d(j)+1.0/(x(j)-x(jj))
        end do
        do jj=j+1,n
          d(j)=d(j)+1.0/(x(j)-x(jj))
        end do
      end do

      do j=0,n
        h(j)=L(j)**2 *(1.0-2.0*(a-x(j))*d(j))
        g(j)=(a-x(j))*L(j)**2
      end do

C     �A���S���Y���菇�R�D
      HH=0.0
      do j=0,n
        HH=HH +f(j)*h(j) +ff(j)*g(j)
      end do

C     ���ʂ̏o��
      write(*,901) "H(0.15)=", HH
      write(*,901) "sin(0.15)=", sin(0.15)
      write(*,901) "error=", sin(0.15)-HH

  901 format(a10,f15.8)
      end
