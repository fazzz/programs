C**********************************************************************
C 1.3�@�ۂߌ덷�Ƒł��؂�덷
C**********************************************************************
      implicit none
      integer :: n,j,k
      real :: S,fact

C n���őł��؂�
      n=10

C �A���S���Y���菇�P�D
      S=1.0

C �A���S���Y���菇�Q�D
      do j=n,1,-1  !�a���E����in,�c,2,1�̏��Ɂj�v�Z����
C       j�̊K��̌v�Z
        fact=1.0
        do k=2,j
          fact=fact*real(k)
        end do

C       S=S+Sj�̌v�Z
        S=S+1.0/fact
      end do

C ���ʂ̏o��
      write(*,900) "calculated:", S, "exact:", exp(1.0)
  900 format(A13, f10.6, A13, f10.6)
      end
