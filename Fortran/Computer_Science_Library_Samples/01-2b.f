C**********************************************************************
C 1.2@‘Q‰»®
C**********************************************************************
      implicit none
      integer :: n,NN
      real :: EPS,x,xx,y,yy,h,gosax,gosay

      NN=700 !Å‘å‚ÅNN‰ñŒJ‚è•Ô‚·
      EPS=1.0E-6 !‘Š‘ÎŒë·‚ª‚±‚ê‚æ‚è¬‚³‚­‚È‚Á‚½‚çû‘©‚µ‚½‚Æ‚·‚é

  900 format(a5,a10,a10)
  901 format(f5.2,f10.6,f10.6) !Œ‹‰Ê‚Ìo—Í‚Ì€”õ
      write(*,900) "nh", "x", "y"

      x=0.0
      y=1.0
      h=0.01
      do n=1,NN
C       ‘Q‰»®‚ÌŒvZ
        xx=x+h*y
        yy=y-h*x
        write(*,901) n*h,x,y

C       û‘©”»’è
        gosax=abs(xx-x)/abs(x)
        gosay=abs(yy-y)/abs(y)
        if(gosax.lt.EPS .and. gosay.lt.EPS) then
          write(*,*) "Converged!"
          exit
        end if

C       Ÿ‚ÌŠÔ‚Ì€”õ
        x=xx
        y=yy
      end do
      end
