C**********************************************************************
C 4.2�@�K�E�X�̏����@�i�Q�j��ޑ��
C**********************************************************************
      implicit none
      integer,parameter :: NM=3
      real :: a(NM,NM),b(NM),x(NM)
      integer :: j,n,k
      real :: sum

      a(1,1)=1.0;   a(1,2)=-4.0;   a(1,3)= 3.0;   b(1)=-1.0
      a(2,1)=0.0;   a(2,2)=-1.0;   a(2,3)=-1.0;   b(2)= 3.0
      a(3,1)=0.0;   a(3,2)= 0.0;   a(3,3)=-5.0;   b(3)=10.0
      n=3

      do j=n,1,-1

        if(j+1.gt.n) then
C         k�i�܂�j+1�ȉ��j>n�̏ꍇ�͑��a�̌v�Z�����Ȃ�
          sum=0.0
        else
C         ���a�̌v�Z
          sum=0.0
          do k=j+1,n
            sum=sum+a(j,k)*x(k)
          end do
        end if

        x(j)=1.0/a(j,j)*(b(j)-sum)

C       ���ʂ̏o��
        write(*,901) "x(", j, ")=", x(j)
      end do
  901 format(a2,I2,a2,f10.6)
      end
