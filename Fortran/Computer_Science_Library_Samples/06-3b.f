C**********************************************************************
C 6.3�@����ȑ�^�A��1���������iSOR�@�j
C**********************************************************************
      implicit none
      integer,parameter :: NI=100,NJ=100
      integer :: i,j,II,JJ,m,mmax
      real :: u(NI,NJ)
      real :: gosa,EPS,us,uu,omg

C ���̐ݒ�
      II=100
      JJ=100
      EPS=1.0E-3
      omg=1.2 !0�`2�̊Ԃ̒l�D1���ƃK�E�X�E�U�C�f���@�Ɉ�v

C �����l
      do j=1,JJ
        do i=1,II
          u(i,j)=0.0
        end do
      end do

C ���͂̉��x
      do j=1,JJ
        u( 1,j)=0.5
        u(II,j)=0.0
      end do
      do i=1,II
        u(i, 1)=1.0
        u(i,JJ)=0.0
      end do

C ������
      mmax=20000
      do m=1,mmax !�����̂��߂�do��
        gosa=0.0

        do j=2,JJ-1
          do i=2,II-1
            us=u(i,j) !�˒i�K�̒l���L���i�덷����Ɏg���j
            uu=(u(i+1,j)+u(i,j+1)+u(i-1,j)+u(i,j-1))/4.0
            gosa=gosa+abs(uu-us)
            u(i,j)=(1.0-omg)*u(i,j)+omg*uu
          end do
        end do

C       �����������ǂ����̔���
        if(gosa.lt.EPS) then
          write(*,*) "Converged!"
          write(*,*) " "
          exit !���[�v���甲����
        else
          if(mod(m,1000).eq.0) write(*,901) "m=", m, "gosa=", gosa
        end if
      end do

C     mmax��J��Ԃ��Ă��������Ȃ������ꍇ
      if(m.gt.mmax) then
        write(*,*) " "
        write(*,*) "It didn't converged."
      end if

C     ���ʂ̏o��
      open(17,file='06-3b.txt', status='replace')
      do j=1,JJ
        do i=1,II
          write(17,902) 0.01*real(i), 0.01*real(j), u(i,j)
        end do
      end do
      close(17)

  901 format(a2,i6,5x,a5,f10.6)
  902 format(2f5.2,f10.6)
      end

