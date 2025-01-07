C**********************************************************************
C 7.3@ƒ„ƒRƒr–@
C**********************************************************************
      implicit none
      integer,parameter :: NM=3
      integer :: i,j,k,kmax,n,p,q
      real :: A(NM,NM),U(NM,NM),atmp(NM,NM),UT(NM,NM),AA(NM,NM)
      real :: EPS,PAI,amax,theta,sum1,sum2,lambda,lambda2

      EPS=1.0E-6
      PAI=4.0*atan(1.0)

C ƒAƒ‹ƒSƒŠƒYƒ€è‡‚PD
      n=3
      a(1,1)=1.0;    a(1,2)=2.0;    a(1,3)=1.0
      a(2,1)=a(1,2); a(2,2)=1.0;    a(2,3)=0.0
      a(3,1)=a(1,3); a(3,2)=a(2,3); a(3,3)=1.0

C ƒAƒ‹ƒSƒŠƒYƒ€è‡‚QD
      kmax=10
      do k=1,kmax

C       ƒAƒ‹ƒSƒŠƒYƒ€è‡‚QD‚P
        amax=0.0
        do i=1,n-1
          do j=i+1,n
            if(abs(a(i,j)).gt.amax) then
              amax=a(i,j)
              p=i
              q=j
            end if
          end do
        end do

C       ƒAƒ‹ƒSƒŠƒYƒ€è‡‚QD‚Q
        if(a(p,p).ne.a(q,q)) then
          theta=0.5*atan(2.0*a(p,q)/(a(q,q)-a(p,p)))
          if(theta.lt.-PAI/2.0 .or. theta.gt.PAI/2.0) then
            write(*,*) "fear of error: theta=", theta
          end if
        else
          theta=sign(PAI/4.0,a(p,q))
          !sign(a1,a2)‚ÍCa2>=0‚È‚ç|a1|Ca2<0‚È‚ç -|a1|
        end if

C       ƒAƒ‹ƒSƒŠƒYƒ€è‡‚QD‚R
        do j=1,n
          do i=1,n
            if(i.eq.j) then                !‘ÎŠp—v‘f‚Ìê‡
              if(i.eq.p .or. i.eq.q) then  !U(p,p)‚ÆU(q,q)‚ÍcosƒÆ
                U(i,j)=cos(theta)
              else                         !‚»‚êˆÈŠO‚Ì‘ÎŠp—v‘f‚Í1.0
                U(i,j)=1.0
              end if
            else                           !‘ÎŠp—v‘f‚Å‚È‚¢ê‡
              if(i.eq.p .and. j.eq.q) then !U(p,q)‚ÍsinƒÆ
                U(i,j)=sin(theta)
              else if(i.eq.q .and. j.eq.p) then !U(q,p)‚Í-sinƒÆ
                U(i,j)=-sin(theta)
              else                         !‚»‚êˆÈŠO‚Í0.0
                U(i,j)=0.0
              end if
            end if
          end do
        end do

C      iQljŠî–{‰ñ“]s—ñ‚ğo—Í
        write(*,*) "Basic rotation matrix", k
        do i=1,n
          write(*,901) U(i,1),U(i,2),U(i,3)
        end do
        write(*,*) " "

C       ƒAƒ‹ƒSƒŠƒYƒ€è‡‚QD‚S
        call transposition(n,U,UT) !U‚Ì“]’us—ñUT‚Ìì¬
!        call multiplication(n,A,U,Atmp) !A*U‚ÌŠ|‚¯Z„Atmp‚ÉŠi”[
!        call multiplication(n,UT,Atmp,A) !UT*Atmp‚ÌŠ|‚¯Z„V‚½‚ÈA‚ÉŠi”[
        call multiplication(n,UT,A,Atmp) !A*U‚ÌŠ|‚¯Z„Atmp‚ÉŠi”[
        call multiplication(n,Atmp,U,A) !UT*Atmp‚ÌŠ|‚¯Z„V‚½‚ÈA‚ÉŠi”[



C      iQljA‚ğo—Í
        write(*,*) "A", k
        do i=1,n
          write(*,901) A(i,1),A(i,2),A(i,3)
        end do
        write(*,*) " "





C       ƒAƒ‹ƒSƒŠƒYƒ€è‡‚RD
        if(amax.lt.EPS) then !”ñ‘ÎŠp—v‘f‚Ìâ‘Î’lÅ‘å‚Ì‚à‚Ìiamaxj
          write(*,*) "Converged!"
          write(*,*) " "
          exit !ƒ‹[ƒv‚©‚ç”²‚¯‚é
        end if

      end do

C     Œ‹‰Ê‚Ìo—Í
      do i=1,n
        write(*,902) "Lambda(", i, ")=", A(i,i)
      end do

  901 format(3f10.6)
  902 format(a7,I2,a2,f10.6)
      end

C**********************************************************************
C “]’us—ñ‚ğì¬‚·‚éƒTƒuƒ‹[ƒ`ƒ“iU‚Ì“]’us—ñ‚ğUT‚ÉŠi”[‚µ‚Ä•Ô‚·j
C**********************************************************************
      subroutine transposition(n,U,UT)
      implicit none
      integer :: n,i,j
      real :: U(n,n),UT(n,n)

      do j=1,n
        do i=1,n
          UT(j,i)=U(i,j)
        end do
      end do
      end subroutine
C**********************************************************************
C s—ñ‚ÌŠ|‚¯Z‚ğ‚·‚éƒTƒuƒ‹[ƒ`ƒ“iA*B‚ğC‚ÉŠi”[‚µ‚Ä•Ô‚·j
C**********************************************************************
      subroutine multiplication(n,A,B,C)
      implicit none
      integer :: n,i,j,k
      real :: A(n,n),B(n,n),C(n,n)

      do j=1,n
        do i=1,n
          C(i,j)=0.0
          do k=1,n
            C(i,j)=C(i,j)+A(i,k)*B(k,j)
          end do
        end do
      end do
      end subroutine