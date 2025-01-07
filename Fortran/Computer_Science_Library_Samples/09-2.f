C**********************************************************************
C 9.2�@�ŏ����@
C**********************************************************************
      implicit none
      integer,parameter :: NM=10 !�ő�f�[�^�g��
      integer,parameter :: MM=5  !m-1�����̍ő�
      integer :: n,m,j,k
      real :: x(NM),y(NM),S(0:MM*2),T(0:MM),a(0:MM)
      real :: xx,yy,dx

C �����f�[�^�̎w��
      x(1)=0.0
      y(1)=1.01
      x(2)=0.2
      y(2)=1.98
      x(3)=0.4
      y(3)=3.33
      x(4)=0.6
      y(4)=4.02
      x(5)=0.8
      y(5)=4.95
      n=5 !�f�[�^��
      m=2 !m-1�����ŋߎ�����

C �ʂ̗�
      x(1)=0.1
      y(1)=2.4824
      x(2)=0.2
      y(2)=1.9975
      x(3)=0.3
      y(3)=1.6662
      x(4)=0.4
      y(4)=1.3775
      x(5)=0.5
      y(5)=1.0933
      x(6)=0.6
      y(6)=0.7304
      x(7)=0.7
      y(7)=0.4344
      x(8)=0.8
      y(8)=0.2981
      x(9)=0.9
      y(9)=-0.0017
      x(10)=1.0
      y(10)=-0.0026
      n=10 !�f�[�^��
      m=3  !m-1�����ŋߎ�����

C Sj��Tj�̌v�Z
      do j=0,2*(m-1)
        S(j)=0.0
        do k=1,n
          S(j)=S(j)+x(k)**j
        end do
      end do

      do j=0,m-1
        T(j)=0.0
        do k=1,n
          T(j)=T(j)+y(k)*x(k)**j
        end do
      end do

C �A��m+1��1��������������
      call GAUSS(S,T,m,a) !�z��a�Ɍ��ʂ��󂯎��

C ���ʂ̏o�́i����ꂽm-1������21�_���̃O���t�Ƃ��ďo�́j
      dx=x(n)/20.0
      do k=1,21
        xx=x(1)+dx*real(k-1)

        !m=2(1����)�Ȃ�y = a(1)*x + a(0)
        !m=3(2����)�Ȃ�y = a(2)*x*x + a(1)*x + a(0)
        yy=0.0
        do j=m-1,0,-1
          yy=yy+a(j)*xx**j
        end do
        write(*,'(f4.2,f10.6)') xx, yy
      end do

C [�Q�l]����ꂽ��������\������
      write(*,*) " "
      write(*,'(a18)',advance='no') "[for reference] y="
      do j=m-1,0,-1
        if(j.eq.0) then
          write(*,'(a2,f8.4)',advance='no') "+",a(j)
        else if(j.eq.1) then
          write(*,'(a2,f8.4,a2)',advance='no') "+",a(j),"x"
        else
          write(*,'(a2,f8.4,a3,I1)',advance='no') "+",a(j),"x^",j
        end if
      end do
      write(*,*) " "
      end

C**********************************************************************
C �A���������������T�u���[�`���i4.1�߃K�E�X�̏����@�𗘗p�j
C**********************************************************************
      subroutine GAUSS(S,T,m,x) !�z��x��Ԃ�
      implicit none
      integer :: m,i,j,k,l
      real :: S(0:m*2),T(0:m)
      real :: a(m,m),aa(m,m),c(m,m),b(m),bb(m),x(0:m)
      real :: sum

      !�W���s��̎w��
      do j=1,m
        do i=1,m
          a(i,j)=S(i+j-2) !S��T�͔z��̃C���f�b�N�X��0�X�^�[�g
          b(i)=T(i-1)
        end do
      end do

      do l=1,m-1
        do j=l+1,m
          c(j,l)=a(j,l)/a(l,l)
          do k=l,m
            aa(j,k)=a(j,k)-c(j,l)*a(l,k)
          end do
          bb(j)=b(j)-c(j,l)*b(l)
        end do

        do j=l+1,m
          do k=l,m
            a(j,k)=aa(j,k)
          end do
          b(j)=bb(j)
        end do
      end do

      do j=m,1,-1
        if(j+1.gt.m) then
          sum=0.0
        else
          sum=0.0
          do k=j+1,m
            sum=sum+a(j,k)*x(k)
          end do
        end if
        x(j)=1.0/a(j,j)*(b(j)-sum)
      end do

      do j=0,m
        x(j)=x(j+1) !�C���f�b�N�X��0�X�^�[�g�ɂ��ĕԂ�
      end do

      end subroutine
