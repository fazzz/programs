C**********************************************************************
C 14.2�@���̕��@��1�����M�`��������������
C**********************************************************************
      implicit none
      integer,parameter :: NM=10
      real :: dt,dx,c
      real :: x(0:NM),u(0:NM),f(0:NM),uu(0:NM)
      integer :: N,M,j,k

C �A���S���Y���菇�P�D
      N=10  !��Ԃ�10����
      M=100 !t=dt*M=0.2�܂Ōv�Z
      dt=0.002
      dx=1.0/real(N)
      c=dt/(dx)**2 !�v�Z�Ɏg���W�����܂Ƃ߂Ă���

C ���U�_�̒�`
      do j=0,N
        x(j)=real(j)*dx
      end do

C �֐�f(x)�̒�`�i�z��Ƃ��ĕێ��j
      do j=0,N
        if(x(j).le.0.5) then
          f(j)=2.0*x(j)
        else
          f(j)=2.0*(1.0-x(j))
        end if
      end do

C �A���S���Y���菇�Q�D
      do j=0,N
        u(j)=f(j)
      end do      

C ���ʂ̏o�͂̏���
      write(*,*) "    t|  x(0)  x(1) ... x(N)"
  901 format(f6.3)

C �A���S���Y���菇�R�D
      do k=1,M
C       �R�|�P�D���E����
        u(0)=0.0
        u(N)=0.0
C       �R�|�Q�D�v�Z
        do j=1,N-1
          uu(j)=u(j)+c*(u(j-1)-2.0*u(j)+u(j+1))
        end do

C       ���ʂ̏o��
        write(*,901,advance='no') dt*real(k) !t
        write(*,'(a)',advance='no') "| "
        do j=0,N
          write(*,901,advance='no') u(j) !���s�Ȃ�
        end do
        write(*,*) " " !1�s�o�͏I���������s

C       uu(j) ���Ȃ킿 n+1 �̎��̒l�����x�� u(j) �ɃR�s�[�D
C       ���̃X�e�b�v�ł܂� u(j) �̒l�i�����_�ł�n+1�j��p����
C       uu(j)�̒l�i�����_�ł�n+2�j�̒l���v�Z����D
        do j=0,N
          u(j)=uu(j)
        end do
      end do
      end
