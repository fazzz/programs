C**********************************************************************
C 4.4�@�|���o���@
C**********************************************************************
      implicit none
      integer,parameter :: NM=3
      real :: a(NM,NM),aa(NM,NM),m(NM,NM),b(NM),bb(NM),atmp(NM),x(NM)
      integer :: j,k,l,n,kk,i
      real :: amax,btmp,sum

      a(1,1)=1.0;   a(1,2)=-4.0;   a(1,3)=3.0;   b(1)=-1.0
      a(2,1)=1.0;   a(2,2)=-5.0;   a(2,3)=2.0;   b(2)= 2.0
      a(3,1)=1.0;   a(3,2)=-1.0;   a(3,3)=1.0;   b(3)= 0.0
      n=3

      do l=1,n
        write(*,*) "l=",l

C**************************
C       �����s�{�b�g�I��
C**************************
C       ��Βl���ő��a(l,l)����T��
        amax=0.0
        do kk=l,n
          if(abs(a(kk,l)).gt.amax) then
            amax=a(kk,l)
            i=kk
          end if
        end do
        write(*,902) " pibot is", amax, " @i=", i

C       ���̓���ւ�
        if(l.ne.i) then
          do kk=l,n
            atmp(kk)=a(l,kk)
          end do
          btmp=b(l)
          do kk=l,n
            a(l,kk)=a(i,kk)
          end do
          b(l)=b(i)
          do kk=l,n
            a(i,kk)=atmp(kk)
          end do
          b(i)=btmp
          write(*,*) " swap:"
          do kk=1,n
            write(*,901) a(kk,1),a(kk,2),a(kk,3),b(kk)
          end do
          write(*,*) " "
         end if

C**************************
C       �|���o���@
C**************************
        do j=1,n
          if(j.ne.l) then

C           �A���S���Y���菇(i)
            m(j,l)=a(j,l)/a(l,l)

C           �A���S���Y���菇(ii)
            do k=l,n
              aa(j,k)=a(j,k)-m(j,l)*a(l,k)
            end do

C           �A���S���Y���菇(iii)
            bb(j)=b(j)-m(j,l)*b(l)
          end if
        end do

        do j=1,n
          if(j.ne.l) then
            do k=l,n
              a(j,k)=aa(j,k)
            end do
            b(j)=bb(j)
          end if
        end do

C       ���ʂ̏o��
        do j=1,n
          write(*,901) a(j,1),a(j,2),a(j,3),b(j)
        end do
        write(*,*) " "
      end do
  901 format(4f10.6)
  902 format(a10,f10.6,a,I2)

C     �A���S���Y���菇(iv)
      do l=1,n
        x(l)=b(l)/a(l,l)
C       ���ʂ̏o��
        write(*,903) "x(", l, ")=", x(l)
      end do
  903 format(a2,I2,a2,f10.6)

      end
