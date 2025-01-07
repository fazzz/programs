C**********************************************************************
C 9.1�@�X�v���C�����
C**********************************************************************
      implicit none
      integer,parameter :: NM=10
      integer,parameter :: NA=20
      integer :: n,k,i,j
      real :: x(0:NM),f(0:NM),xx(0:NA)
      real :: a(0:NM),b(0:NM),c(0:NM),d(0:NM)
      real :: S2(0:NM)
      real :: SS,Hk

      write(*,902) "xx", "S(xx)"
  902 format(a3,a10)

C �A���S���Y���菇�P�D
      !�l���������Ă��闣�U�_����
      n=8  !��Ԑ�
      do i=0,n
        x(i)=1.0+0.25*real(i)
        f(i)=1.0/(1.0+10.0*x(i)**8)
      end do

      !1����3�܂�0.1���݂�21�_���́uxx�v�ɑ΂��āuS(xx)�v���v�Z����
      do i=0,NA
        xx(i)=1.0+0.1*real(i)
      end do

C �A���S���Y���菇�Q�D
      !S��2�K�����ɑ΂���A��1���������𗧂Ă�
      !3���������ɂȂ�̂Ńg�[�}�X�@�𗘗p
      do k=1,n-1
        a(k)=x(k)-x(k-1)
        b(k)=2.0*(x(k+1)-x(k-1))
        c(k)=x(k+1)-x(k)
        d(k)=6.0*( (f(k+1)-f(k))/(x(k+1)-x(k))
     &            -(f(k)-f(k-1))/(x(k)-x(k-1)))
      end do
      a(1)=0.0
      c(n-1)=0.0

C �A���S���Y���菇�R�D
      S2(0)=0.0
      S2(n)=0.0

C �A���S���Y���菇�S�D
      call THOMAS(1,n-1,a,b,c,d,S2)

C �exx(i)�ɑ΂��āA�ǂ̏���Ԃɂ���̂�
C ���Ȃ킿�A�ǂ̕�Ԏ�S(k)���g���ׂ����T��
      do i=0,NA
        if(xx(i).lt.x(0)) then
          k=0
        else if(x(n).le.xx(i)) then
          k=n-1
        else
          do j=0,n-1
          if(x(j).le.xx(i) .and. xx(i).lt.x(j+1)) then
            k=j
            exit
          end if
          end do
        end if

C �A���S���Y���菇�T�D
        Hk=x(k+1)-x(k)
        SS= (S2(k+1)*(xx(i)-x(k))**3 - S2(k)*(xx(i)-x(k+1))**3)/(6.0*Hk)
     &     +(f(k+1)/Hk - Hk*S2(k+1)/6.0)*(xx(i)-x(k))
     &     -(f(k)/Hk - Hk*S2(k)/6.0)*(xx(i)-x(k+1))
        write(*,901) xx(i),SS
      end do
  901 format(f3.1,f10.6)

      end

C**********************************************************************
C �g�[�}�X�@���v�Z���邽�߂̃T�u���[�`��
C**********************************************************************
      subroutine THOMAS(IL,IU,a,b,c,d,y)
      implicit none
      integer,parameter :: NM=10
      integer :: IL,IU,j
      real :: a(0:NM),b(0:NM),c(0:NM),d(0:NM)
      real :: g(0:NM),s(0:NM)
      real :: x(0:NM),y(0:NM)

      g(IL)=b(IL)
      s(IL)=d(IL)
      do j=IL+1,IU
        g(j)=b(j)-a(j)*c(j-1)/g(j-1)
        s(j)=d(j)-a(j)*s(j-1)/g(j-1)
      end do
      y(IU)=s(IU)/g(IU)
      do j=IU-1,IL,-1
        y(j)=(s(j)-c(j)*y(j+1))/g(j)
      end do
      end subroutine
