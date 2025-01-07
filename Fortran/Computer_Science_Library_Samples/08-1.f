C**********************************************************************
C 8.1�@���O�����W�����
C**********************************************************************
      implicit none
      integer,parameter :: NM=10
      integer,parameter :: NA=20
      integer :: n,j,jj,k
      real :: x(0:NM),f(0:NM),l(0:NM),a(0:NA)
      real :: L_bunshi1,L_bunshi2,L_bunbo1,L_bunbo2,P

      write(*,902) "a", "P(a)"
  902 format(a3,a10)

C �A���S���Y���菇�P�D
      !�l���������Ă��闣�U�_�����̎w��
      n=3  !3�����̏ꍇ�i4�_�j
      x(0)=0.0
      f(0)=1.00000
      x(1)=1.0
      f(1)=0.09091
      x(2)=2.0
      f(2)=0.00039
      x(3)=3.0
      f(3)=0.00002

      n=6  !6�����̏ꍇ�i7�_�j
      x(0)=0.0
      f(0)=1.00000
      x(1)=0.5
      f(1)=0.96241
      x(2)=1.0
      f(2)=0.09091
      x(3)=1.5
      f(3)=0.00389
      x(4)=2.0
      f(4)=0.00039
      x(5)=2.5
      f(5)=0.00007
      x(6)=3.0
      f(6)=0.00002

      !1����3�܂�0.1���݂�21�_���́ua�v�ɑ΂��āuP(a)�v���v�Z����
      do k=0,NA
        a(k)=1.0+0.1*real(k)

C       �A���S���Y���菇�Q�D
        do j=0,n
          L_bunshi1=1.0
          L_bunshi2=1.0
          L_bunbo1 =1.0
          L_bunbo2 =1.0
          do jj=0,j-1
            L_bunshi1=L_bunshi1*(a(k)-x(jj))   !(a-X0)(a-X1)�c(a-Xj-1)
            L_bunbo1 =L_bunbo1 *(x(j)-x(jj))   !(Xj-X0)(Xj-X1)�c(Xj-Xj-1)
          end do
          do jj=j+1,n
            L_bunshi2=L_bunshi2*(a(k)-x(jj))   !(a-Xj+1)�c(a-Xn)
            L_bunbo2 =L_bunbo2 *(x(j)-x(jj))   !(Xj-Xj+1)�c(Xj-Xn)
          end do
          L(j)=(L_bunshi1*L_bunshi2)/(L_bunbo1*L_bunbo2)
        end do

C       �A���S���Y���菇�R�D
        P=0.0
        do j=0,n
          P=P+f(j)*l(j)
        end do
        write(*,901) a(k), P

      end do !do k �ɑ΂���end do
  901 format(f3.1,f10.6)

      end
