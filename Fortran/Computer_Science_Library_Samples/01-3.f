C**********************************************************************
C 1.3　丸め誤差と打ち切り誤差
C**********************************************************************
      implicit none
      integer :: n,j,k
      real :: S,fact

C n項で打ち切る
      n=10

C アルゴリズム手順１．
      S=1.0

C アルゴリズム手順２．
      do j=n,1,-1  !和を右から（n,…,2,1の順に）計算する
C       jの階乗の計算
        fact=1.0
        do k=2,j
          fact=fact*real(k)
        end do

C       S=S+Sjの計算
        S=S+1.0/fact
      end do

C 結果の出力
      write(*,900) "calculated:", S, "exact:", exp(1.0)
  900 format(A13, f10.6, A13, f10.6)
      end
