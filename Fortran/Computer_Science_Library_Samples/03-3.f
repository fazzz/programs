C**********************************************************************
C 3.3　ベアストウ法
C**********************************************************************
      implicit none
      integer,parameter :: NM=10
      real :: a(0:NM),b(-2:NM),c(-2:NM)
      real :: EPS,u,v,du,dv,uu,vv
      integer :: k,n,i,j,JJ,NN

C     係数aの指定
!      a(0)=1.0;  a(1)=3.0; a(2)=-15.0; a(3)=-24.0; a(4)=15.0; a(5)=63.0
!      a(6)=71.0; a(7)=30.0
!      NN=7

      a(0)=1.0; a(1)=0.0; a(2)=-5.0; a(3)=-9.0; a(4)=-8.0; a(5)=-3.0
      NN=5

!      a(0)=1.0; a(1)=3.0; a(2)=4.0; a(3)=3.0; a(4)=1.0
!      NN=4

!      a(0)=1.0; a(1)=-9.0; a(2)=23.0; a(3)=-15.0
!      a(0)=1.0; a(1)=0.0; a(2)=0.0; a(3)=-1.0
!      a(0)=1.0; a(1)=-3.0; a(2)=3.0; a(3)=-1.0
!      NN=3

!      a(0)=1.0; a(1)=1.0; a(2)=-6.0
!      NN=2

      EPS=1.0E-6
      JJ=50
      n=NN

      do i=0,NN

C       初期値
        uu=1.0 !uの近似値(.ne.0)
        vv=1.0 !vの近似値(.ne.0)

        do j=1,JJ  !連立方程式の収束のための反復
C         アルゴリズム手順１．
          do k=0,n
            b(-2)=0.0
            b(-1)=0.0
            b(k)=a(k)+uu*b(k-1)+vv*b(k-2)
          end do

C         アルゴリズム手順２．
          do k=0,n-1
            c(-2)=0.0
            c(-1)=0.0
            c(k)=+b(k)+uu*c(k-1)+vv*c(k-2)
          end do

C         アルゴリズム手順３．（3.2節を利用）
          call simultaneous(b(n-1),b(n),c(n-2),c(n-3),c(n-1),c(n-2),
     &                      du,dv) !du,dvに計算結果を記憶
          u=uu+du
          v=vv+dv

C         アルゴリズム手順４．
          if(abs(du).lt.EPS .and. abs(dv).lt.EPS) then
            exit !終了
          else
            uu=u
            vv=v
          end if
        end do

C       アルゴリズム手順５．
        call quadric(1.0,-u,-v)

C       アルゴリズム手順６．
        n=n-2
        if(n.eq.0) then
          exit
        else if(n.eq.1) then
          write(*,903) "x=", -b(1)/b(0)
          exit
        else if(n.eq.2) then
          call quadric(b(0),b(1),b(2))
          exit
        else
          do k=0,n
            a(k)=b(k)
          end do
        end if
      end do
  903 format(a2,f10.6)

      end

C**********************************************************************
C 二次方程式の解を求めるサブルーチン（1.4節を利用）
C**********************************************************************
      subroutine quadric(a,b,c)
      implicit none
      real :: a,b,c,D,x1,RE,IM

  900 format(a2,f10.6,a3,f10.6,a1)
  901 format(a2,f10.6,a)
  902 format(a2,f10.6)

      D=b**2.0-4.0*a*c
      if(D.lt.0.0) then
        RE=-b/(2.0*a)
        IM=sqrt(abs(D))/(2.0*a)
        write(*,900) "x=",RE,"+",IM,"i"
        write(*,900) "x=",RE,"-",IM,"i"
      else if(D.eq.0.0) then
        write(*,901) "x=", -b/(2.0*a), "(multiple root)"
      else
        if(b.ge.0.0) then
          x1=(-b-sqrt(D))/(2.0*a)
          write(*,902) "x=", x1
          write(*,902) "x=", c/(a*x1)
        else
          x1=(-b+sqrt(D))/(2.0*a)
          write(*,902) "x=", x1
          write(*,902) "x=", c/(a*x1)
        end if
      end if
      end subroutine

C**********************************************************************
C 連立方程式を解くサブルーチン（3.2節を利用）
C**********************************************************************
      subroutine simultaneous(f,g,fx,fy,gx,gy,dx,dy)
      implicit none
      real :: f,g,fx,fy,gx,gy,dx,dy,aa

      aa=fy*gx-fx*gy
      if(aa.ne.0.0) then
        dx=(fy*(-g)-gy*(-f))/aa
        dy=(gx*(-f)-fx*(-g))/aa
      else
        write(*,*) "aa=0.0 (fear of error)"
        dx=1.0 !ダミーの値
        dy=1.0
      end if
      end subroutine
