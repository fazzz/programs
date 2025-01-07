C**********************************************************************
C 4.1　ガウスの消去法（１）前進消去
C**********************************************************************
      implicit none
      integer,parameter :: NM=3
      real :: a(NM,NM),aa(NM,NM),m(NM,NM),b(NM),bb(NM)
      integer :: j,k,l,n

      a(1,1)=1.0;   a(1,2)=-4.0;   a(1,3)=3.0;   b(1)=-1.0
      a(2,1)=1.0;   a(2,2)=-5.0;   a(2,3)=2.0;   b(2)= 2.0
      a(3,1)=1.0;   a(3,2)=-1.0;   a(3,3)=1.0;   b(3)= 0.0
      n=3

      do l=1,n-1
        write(*,*) "l=",l

        do j=l+1,n
C         アルゴリズム手順(i)
          m(j,l)=a(j,l)/a(l,l)

C         アルゴリズム手順(ii)
          do k=l,n !本質的にはdo k=l+1,nでよい
            aa(j,k)=a(j,k)-m(j,l)*a(l,k)
          end do

C         アルゴリズム手順(iii)
          bb(j)=b(j)-m(j,l)*b(l)
        end do

        do j=l+1,n
          do k=l,n !本質的にはdo k=l+1,nでよい
            a(j,k)=aa(j,k)
          end do
          b(j)=bb(j)
        end do

C       結果の出力
        do j=1,n
          write(*,901) a(j,1),a(j,2),a(j,3),b(j)
        end do
        write(*,*) " "
      end do
  901 format(4f10.6)
      end
