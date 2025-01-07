C**********************************************************************
C 7.2�@�t�x�L��@
C**********************************************************************
      implicit none
      integer,parameter :: NM=3
      integer :: i,j,k,kmax,n
      real :: a(NM,NM),x(NM),y(NM),z(NM),U(NM,NM),L(NM,NM)
      real :: EPS,sum1,sum2,lambda,lambda2,gosa

      EPS=1.0E-6

C �A���S���Y���菇�P�D
      n=3
      a(1,1)=11.0; a(1,2)= 7.0; a(1,3)=-5.0
      a(2,1)= 0.0; a(2,2)=10.0; a(2,3)=-1.0
      a(3,1)= 2.0; a(3,2)= 8.0; a(3,3)= 3.0

      !�K���ł悢���A�[���x�N�g���ɂ���ƌv�Z�ł��Ȃ�
      do i=1,n
        x(i)=1.0
      end do

C �A���S���Y���菇�Q�D
      call LUmatrix(a,U,L,n)

C �A���S���Y���菇�R�D
      kmax=200
      lambda=10000.0 !�ɂ̏����l�i�_�~�[�j
      lambda2=10000.0
      do k=1,kmax
        gosa=0.0

C       (1)LUy=x������
C       �܂�Lz=x������
        z(1)=x(1)
        do j=2,n
          sum1=x(j)
          do i=1,j-1
            sum1=sum1-L(j,i)*z(i)
          end do
          z(j)=sum1
        end do
C       �Â���Uy=z������
        y(n)=z(n)/U(n,n)
        do j=n-1,1,-1
          sum1=z(j)
          do i=j+1,n
            sum1=sum1-U(j,i)*y(i)
          end do
          y(j)=sum1/U(j,j)
        end do

C       ��=(y,x)/(y,y)
        sum1=0.0
        sum2=0.0
        do i=1,n
          sum1=sum1+y(i)*x(i)
          sum2=sum2+y(i)*y(i)
        end do
        lambda=sum1/sum2

C       (2)
        do i=1,n
          x(i)=y(i)/sqrt(sum2)
        end do

C       �A���S���Y���菇�S�D
        gosa=abs(lambda-lambda2)/abs(lambda2)
        lambda2=lambda
        if(gosa.lt.EPS) then
          write(*,*) "Converged!"
          write(*,*) " "
          exit !���[�v���甲����
        end if

C       ���ʂ̏o��
        write(*,901) k,lambda
      end do

C     kmax��J��Ԃ��Ă��������Ȃ������ꍇ
      if(k.gt.kmax) then
        write(*,*) " "
        write(*,*) "It didn't converged."
      end if

  901 format(I3,f10.6)
      end

C**********************************************************************
C LU����������T�u���[�`��
C**********************************************************************
      subroutine LUmatrix(a,U,L,n)
      implicit none
      integer :: i,j,k,n
      real :: aa(n,n),a(n,n),U(n,n),L(n,n)

C     �z��L�̏�����
      do j=1,n
        do i=1,n
          L(i,j)=0.0
        end do
      end do
      do i=1,n
        L(i,i)=1.0
      end do

      do j=1,n
        do i=1,n
          U(i,j)=a(i,j)
        end do
      end do

      do j=1,n-1
        do i=j+1,n
          L(i,j)=U(i,j)/U(j,j)
          do k=j,n
            aa(i,k)=U(i,k)-L(i,j)*U(j,k)
          end do
        end do
        do i=j+1,n
          do k=j,n
            U(i,k)=aa(i,k)
          end do
        end do
      end do

      write(*,*) " Upper triangular matrix"
      do i=1,n
        write(*,901) U(i,1),U(i,2),U(i,3)
      end do
      write(*,*) " "
      write(*,*) " Lower triangular matrix"
      do i=1,n
        write(*,901) L(i,1),L(i,2),L(i,3)
      end do
      write(*,*) " "
  901 format(4f10.6)
      end subroutine
