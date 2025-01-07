C**********************************************************************
C 11.4�@���U�t�[���G�ϊ�
C**********************************************************************
      implicit none
      integer,parameter :: NM=20 !n�ȏ�ɂ��邱��
      integer :: n,k,j
      real :: PAI,T,c,phi_real,phi_imag,dx,u
      real :: w_real(0:NM,0:NM),w_imag(0:NM,0:NM),x(0:NM)
      real :: f

      PAI=4.0*atan(1.0)
      n=16 !�f�[�^��
      dx=0.25 !�T���v�����O�Ԋu
      T=dx*real(n)

      write(*,900) "k", "real part", "imaginary part"

C �A���S���Y���菇�P�D
      do k=0,n-1

C       �A���S���Y���菇�P�D�P
        do j=0,n-1
          c=2.0*PAI*real(k)*real(j)/real(n)
          w_real(k,j)=cos(c)
          w_imag(k,j)=-sin(c)
        end do

C       �A���S���Y���菇�P�D�Q
        phi_real=0.0
        phi_imag=0.0
        do j=0,n-1
          x(j)=T*real(j)/real(n)
          phi_real=phi_real+w_real(k,j)*f(x(j))
          phi_imag=phi_imag+w_imag(k,j)*f(x(j))
        end do

        write(*,901) k,phi_real*dx,phi_imag*dx
      end do



       write(*,*) " "
      do k=0,n-1
        u=PAI*2.0*real(k)

        write(*,901) k, 1.0/(1+u**2), -u/(1+u**2)
      end do
















  900 format(A3,A20,A20)
  901 format(I3,f20.8,f20.8)
      end

C**********************************************************************
C �֐�f�̒�`
C**********************************************************************
      real function f(t)
      implicit none
      real :: t

      f=1.0/exp(t)
      end function
