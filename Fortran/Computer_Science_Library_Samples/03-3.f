C**********************************************************************
C 3.3�@�x�A�X�g�E�@
C**********************************************************************
      implicit none
      integer,parameter :: NM=10
      real :: a(0:NM),b(-2:NM),c(-2:NM)
      real :: EPS,u,v,du,dv,uu,vv
      integer :: k,n,i,j,JJ,NN

C     �W��a�̎w��
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

C       �����l
        uu=1.0 !u�̋ߎ��l(.ne.0)
        vv=1.0 !v�̋ߎ��l(.ne.0)

        do j=1,JJ  !�A���������̎����̂��߂̔���
C         �A���S���Y���菇�P�D
          do k=0,n
            b(-2)=0.0
            b(-1)=0.0
            b(k)=a(k)+uu*b(k-1)+vv*b(k-2)
          end do

C         �A���S���Y���菇�Q�D
          do k=0,n-1
            c(-2)=0.0
            c(-1)=0.0
            c(k)=+b(k)+uu*c(k-1)+vv*c(k-2)
          end do

C         �A���S���Y���菇�R�D�i3.2�߂𗘗p�j
          call simultaneous(b(n-1),b(n),c(n-2),c(n-3),c(n-1),c(n-2),
     &                      du,dv) !du,dv�Ɍv�Z���ʂ��L��
          u=uu+du
          v=vv+dv

C         �A���S���Y���菇�S�D
          if(abs(du).lt.EPS .and. abs(dv).lt.EPS) then
            exit !�I��
          else
            uu=u
            vv=v
          end if
        end do

C       �A���S���Y���菇�T�D
        call quadric(1.0,-u,-v)

C       �A���S���Y���菇�U�D
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
C �񎟕������̉������߂�T�u���[�`���i1.4�߂𗘗p�j
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
C �A���������������T�u���[�`���i3.2�߂𗘗p�j
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
        dx=1.0 !�_�~�[�̒l
        dy=1.0
      end if
      end subroutine
