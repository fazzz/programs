C**********************************************************************
C 15.3�@�|�A�\���������̍�����@
C**********************************************************************
      implicit none
      integer,parameter :: NM=10
      real :: dx,dy,c
      real :: x(0:NM,0:NM),y(0:NM,0:NM),u(0:NM,0:NM),uu(0:NM,0:NM)
      real :: f(0:NM,0:NM)
      integer :: JJ,J1,KK,K1,j,k,NN,n

C �e�萔�Ȃǂ̒�`
      JJ=10   !x�����̊i�q��
      J1=JJ-1
      KK=10   !y�����̊i�q��
      K1=KK-1
      dx=1.0/real(JJ) !�̈��1.0�~1.0�̐����`
      dy=1.0/real(KK)
      c=1.0/(2.0/(dx)**2+2.0/(dy)**2) !�W���͂܂Ƃ߂Ă���
      NN=50   !�K�E�X�E�U�C�f���@�̔�����

C �֐�f(x)�̒�`�i�z��Ƃ��ĕێ��j
      do k=0,KK
      do j=0,JJ
        f(j,k)=0.0
      end do
      end do

C ��������
      do k=0,KK
      do j=0,JJ
        u(j,k)=0.0
       uu(j,k)=0.0
      end do
      end do
      do j=0,JJ
        u(j,0)=1.0
       uu(j,0)=1.0
      end do

C �K�E�X�E�U�C�f���@�̔���
      do n=1,NN

C       �������̌v�Z
        do k=1,K1
        do j=1,J1
          uu(j,k)=c*( (uu(j-1,k)+u(j+1,k))/(dx)**2
     &               +(uu(j,k-1)+u(j,k+1))/(dy)**2 +f(j,k))
        end do
        end do

C       �v�Z���ʂ̃R�s�[
        do k=1,K1
        do j=1,J1
          u(j,k)=uu(j,k)
        end do
        end do
      end do

C     ���ʂ̏o��
      do k=0,KK
        do j=0,J1
          write(*,901,advance='no') u(j,k) !���s�Ȃ�
        end do
        write(*,901) u(JJ,k) !�Ōゾ�����s
      end do
  901 format(f6.3)

      end
