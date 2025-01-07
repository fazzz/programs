C**********************************************************************
C 4.2　ガウスの消去法（２）部分ピボット選択
C**********************************************************************
      implicit none
      integer,parameter :: NM=3
      real :: a(NM,NM),aa(NM,NM),m(NM,NM),b(NM),bb(NM),atmp(NM)
      integer :: j,k,l,n,kk,i
      real :: amax,btmp

      a(1,1)=1.0;   a(1,2)=-4.0;   a(1,3)=3.0;   b(1)=-1.0
      a(2,1)=1.0;   a(2,2)=-5.0;   a(2,3)=2.0;   b(2)= 2.0
      a(3,1)=1.0;   a(3,2)=-1.0;   a(3,3)=1.0;   b(3)= 0.0
      n=3

      do l=1,n-1
        write(*,*) "l=",l

C       絶対値が最大のa(l,l)候補を探す
        amax=0.0
        do kk=l,n
          if(abs(a(kk,l)).gt.amax) then
            amax=a(kk,l)
            i=kk
          end if
        end do
        write(*,902) " pibot is", amax, " @i=", i

C       式の入れ替え
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
          write(*,*) " swap: "
          do kk=1,n
            write(*,901) a(kk,1),a(kk,2),a(kk,3),b(kk)
          end do
          write(*,*) " "
         end if

C       この後は4.1節と同じ
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
  902 format(a10,f10.6,a,I2)
      end
