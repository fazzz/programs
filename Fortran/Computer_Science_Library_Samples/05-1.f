C**********************************************************************
C 5.1　LU分解法
C**********************************************************************
      implicit none
      integer,parameter :: NM=3
      real :: a(NM,NM),aa(NM,NM),m(NM,NM),b(NM),x(NM),y(NM)
      integer :: j,k,l,n
      real :: sum

      a(1,1)=2.0;   a(1,2)=1.0;   a(1,3)= 1.0;   b(1)=-1.0
      a(2,1)=4.0;   a(2,2)=3.0;   a(2,3)= 4.0;   b(2)= 2.0
      a(3,1)=6.0;   a(3,2)=5.0;   a(3,3)=10.0;   b(3)= 2.0
      n=3

C     配列mの初期化
      do l=1,n
        do k=1,n
          m(k,l)=0.0
        end do
      end do
      do l=1,n
        m(l,l)=1.0
      end do

      do l=1,n
        write(*,*) "l=",l

C       前進消去
        do j=l+1,n
          m(j,l)=a(j,l)/a(l,l)
          do k=l,n
            aa(j,k)=a(j,k)-m(j,l)*a(l,k)
          end do
        end do
        do j=l+1,n
          do k=l,n
            a(j,k)=aa(j,k)
          end do
        end do

        write(*,*) " Upper triangular matrix"
        do j=1,n
          write(*,901) a(j,1),a(j,2),a(j,3)
        end do
        write(*,*) " "
      end do

      write(*,*) " Lower triangular matrix"
      do j=1,n
        write(*,901) m(j,1),m(j,2),m(j,3)
      end do
      write(*,*) " "
  901 format(4f10.6)

C     式(5.9)の1つめの式を解く
      y(1)=b(1)
      do l=2,n
        sum=b(l)
        do j=1,l-1
          sum=sum-m(l,j)*y(j)
        end do
        y(l)=sum
      end do

C     式(5.9)の2つめの式を解く
      x(n)=y(n)/a(n,n)
      do l=n-1,1,-1
        sum=y(l)
        do j=l+1,n
          sum=sum-a(l,j)*x(j)
        end do
        x(l)=sum/a(l,l)
      end do

      do l=1,n
        write(*,903) "x(", l, ")=", x(l)
      end do
  903 format(a2,I2,a2,f10.6)

      end
