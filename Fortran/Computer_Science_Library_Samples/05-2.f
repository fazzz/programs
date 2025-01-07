C**********************************************************************
C 5.2　変形コレスキー法
C**********************************************************************
      implicit none
      integer,parameter :: NM=3
      real :: a(NM,NM),b(NM),x(NM),y(NM),z(NM)
      integer :: j,i,m,n,kk
      real :: sum

C     アルゴリズム手順１．
      a(1,1)=1.0;    a(1,2)=1.0;    a(1,3)= 2.0;   b(1)=-1.0
      a(2,1)=a(1,2); a(2,2)=3.0;    a(2,3)= 4.0;   b(2)= 2.0
      a(3,1)=a(1,3); a(3,2)=a(2,3); a(3,3)=10.0;   b(3)=13.0
      n=3

C     アルゴリズム手順２．
      do j=2,n
        a(1,j)=a(1,j)/a(1,1)
      end do

      do i=2,n
        sum=0.0
        do m=1,i-1
          sum=sum+a(m,i)**2 *a(m,m)
        end do
        a(i,i)=a(i,i)-sum

        do j=i+1,n
          sum=0.0
          do m=1,i-1
            sum=sum+a(m,i)*a(m,m)*a(m,j)
          end do
          a(i,j)=(a(i,j)-sum)/a(i,i)
        end do
      end do

      write(*,*) "D:"
        write(*,902) a(1,1), 0.0,    0.0
        write(*,902) 0.0,    a(2,2), 0.0
        write(*,902) 0.0,    0.0,    a(3,3)
      write(*,*) " "
      write(*,*) "U:"
        write(*,902) 1.0,    a(1,2), a(1,3)
        write(*,902) 0.0,    1.0,    a(2,3)
        write(*,902) 0.0,    0.0,    1.0
      write(*,*) " "

C     アルゴリズム手順３．
      z(1)=b(1)
      do i=2,n
        sum=0.0
        do m=1,i-1
          sum=sum+a(m,i)*z(m)
        end do
        z(i)=b(i)-sum
      end do

C     アルゴリズム手順４．
      x(n)=z(n)/a(n,n)
      do i=n-1,1,-1
        sum=0.0
        do m=i+1,n
          sum=sum+a(i,m)*x(m)
        end do
        x(i)=z(i)/a(i,i)-sum
      end do

      do i=1,n
        write(*,903) "x(", i, ")=", x(i)
      end do

  902 format(3f10.6)
  903 format(a2,I2,a2,f10.6)

      end
