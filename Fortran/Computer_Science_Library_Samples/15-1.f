C**********************************************************************
C 15.1�@�ڗ��������̍�����@
C**********************************************************************
      implicit none
      integer,parameter :: NM=20
      real :: dt,dx,c,r
      real :: x(0:NM),u(0:NM),f(0:NM),uu(0:NM)
      integer :: JJ,NN,j,n

C �e�萔�Ȃǂ̒�`
      JJ=20            !x�����̗��U�_�̐�
      dx=2.0/real(JJ)  !�̈��JJ����
      dt=0.01          !dt�͎w��
      NN=30            !t�����̗��U�_�̐��iNN�X�e�b�v�ڂ܂Ōv�Z�j
      c=2.0            !�ڗ����x
      r=c*dt/dx

C ���U�_�̒�`
      do j=0,JJ
        x(j)=real(j)*dx
      end do

C �֐�f(x)�̒�`�i�z��Ƃ��ĕێ��j
      do j=0,JJ
        f(j)=sin(3.0*x(j))
        if(f(j).lt.0.0) f(j)=0.0
      end do

C ���������i��(15.3)�j
      do j=0,JJ
        u(j)=f(j)
       uu(j)=f(j) 
      end do      
      call out(dt,n,JJ,u) !���������̏o��

C ���Ԃɑ΂���1����NN�܂�1�s���v�Z
      do n=1,NN
C       ��(15.4)�̌v�Z
C       ���Ӂit=n+1�̂Ƃ��̒l�j��z��uu�Ɋi�[����
        do j=1,JJ
          uu(j)=(1-r)*u(j)+r*u(j-1)
        end do

C       �v�Z����uu�̏o��
        call out(dt,n,JJ,uu)

C       uu(j) ���Ȃ킿 n+1 �̎��̒l�� u(j) �ɃR�s�[�D
C       �i���̎��Ԃ̏����j
        do j=0,JJ
          u(j)=uu(j)
        end do
      end do !����Ŏ��Ԃ̌v�Z����s�I��
      end

C**********************************************************************
C ���ʂ��o�͂��邽�߂̃T�u���[�`��
C**********************************************************************
      subroutine out(dt,n,JJ,u)
      real :: dt
      real :: u(JJ)
      integer :: n,j,JJ

      write(*,901,advance='no') dt*real(n) !t
      write(*,'(a)',advance='no') "| "
      do j=0,JJ-1
        write(*,902,advance='no') u(j) !���s�Ȃ�
      end do
      write(*,902) u(JJ) !�Ōゾ�����s

  901 format(f6.2)
  902 format(f6.3)
      end subroutine
