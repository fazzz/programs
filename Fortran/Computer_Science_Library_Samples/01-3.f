C**********************************************************************
C 1.3@Ϋίλ·ΖΕΏΨθλ·
C**********************************************************************
      implicit none
      integer :: n,j,k
      real :: S,fact

C nΕΕΏΨι
      n=10

C ASYθPD
      S=1.0

C ASYθQD
      do j=n,1,-1  !aπE©ηin,c,2,1ΜΙjvZ·ι
C       jΜKζΜvZ
        fact=1.0
        do k=2,j
          fact=fact*real(k)
        end do

C       S=S+SjΜvZ
        S=S+1.0/fact
      end do

C ΚΜoΝ
      write(*,900) "calculated:", S, "exact:", exp(1.0)
  900 format(A13, f10.6, A13, f10.6)
      end
