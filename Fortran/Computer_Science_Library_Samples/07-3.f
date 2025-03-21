C**********************************************************************
C 7.3　ヤコビ法
C**********************************************************************
      implicit none
      integer,parameter :: NM=3
      integer :: i,j,k,kmax,n,p,q
      real :: A(NM,NM),U(NM,NM),atmp(NM,NM),UT(NM,NM),AA(NM,NM)
      real :: EPS,PAI,amax,theta,sum1,sum2,lambda,lambda2

      EPS=1.0E-6
      PAI=4.0*atan(1.0)

C アルゴリズム手順１．
      n=3
      a(1,1)=1.0;    a(1,2)=2.0;    a(1,3)=1.0
      a(2,1)=a(1,2); a(2,2)=1.0;    a(2,3)=0.0
      a(3,1)=a(1,3); a(3,2)=a(2,3); a(3,3)=1.0

C アルゴリズム手順２．
      kmax=10
      do k=1,kmax

C       アルゴリズム手順２．１
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

C       アルゴリズム手順２．２
        if(a(p,p).ne.a(q,q)) then
          theta=0.5*atan(2.0*a(p,q)/(a(q,q)-a(p,p)))
          if(theta.lt.-PAI/2.0 .or. theta.gt.PAI/2.0) then
            write(*,*) "fear of error: theta=", theta
          end if
        else
          theta=sign(PAI/4.0,a(p,q))
          !sign(a1,a2)は，a2>=0なら|a1|，a2<0なら -|a1|
        end if

C       アルゴリズム手順２．３
        do j=1,n
          do i=1,n
            if(i.eq.j) then                !対角要素の場合
              if(i.eq.p .or. i.eq.q) then  !U(p,p)とU(q,q)はcosθ
                U(i,j)=cos(theta)
              else                         !それ以外の対角要素は1.0
                U(i,j)=1.0
              end if
            else                           !対角要素でない場合
              if(i.eq.p .and. j.eq.q) then !U(p,q)はsinθ
                U(i,j)=sin(theta)
              else if(i.eq.q .and. j.eq.p) then !U(q,p)は-sinθ
                U(i,j)=-sin(theta)
              else                         !それ以外は0.0
                U(i,j)=0.0
              end if
            end if
          end do
        end do

C      （参考）基本回転行列を出力
        write(*,*) "Basic rotation matrix", k
        do i=1,n
          write(*,901) U(i,1),U(i,2),U(i,3)
        end do
        write(*,*) " "

C       アルゴリズム手順２．４
        call transposition(n,U,UT) !Uの転置行列UTの作成
!        call multiplication(n,A,U,Atmp) !A*Uの掛け算＞Atmpに格納
!        call multiplication(n,UT,Atmp,A) !UT*Atmpの掛け算＞新たなAに格納
        call multiplication(n,UT,A,Atmp) !A*Uの掛け算＞Atmpに格納
        call multiplication(n,Atmp,U,A) !UT*Atmpの掛け算＞新たなAに格納



C      （参考）Aを出力
        write(*,*) "A", k
        do i=1,n
          write(*,901) A(i,1),A(i,2),A(i,3)
        end do
        write(*,*) " "





C       アルゴリズム手順３．
        if(amax.lt.EPS) then !非対角要素の絶対値最大のもの（amax）
          write(*,*) "Converged!"
          write(*,*) " "
          exit !ループから抜ける
        end if

      end do

C     結果の出力
      do i=1,n
        write(*,902) "Lambda(", i, ")=", A(i,i)
      end do

  901 format(3f10.6)
  902 format(a7,I2,a2,f10.6)
      end

C**********************************************************************
C 転置行列を作成するサブルーチン（Uの転置行列をUTに格納して返す）
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
C 行列の掛け算をするサブルーチン（A*BをCに格納して返す）
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