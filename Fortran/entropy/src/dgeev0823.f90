program bestfit_ty

use xdr, only: xtcfile
use bestfit
        implicit none
        type(xtcfile) :: xtcbf
        type(xtcfile) :: xtc_out
        integer :: i, j, frame, natom
        double precision :: center(3), entropy
        double precision, allocatable, dimension(:) :: stmass
        double precision, allocatable, dimension(:,:) :: ev, varm
        double precision, allocatable, dimension(:,:,:) :: refcoord, cencoord, vec
        character bffilename*50
        real, parameter :: massH = 1.00794          ! atomic weight (hydrogen)
        real, parameter :: massC = 12.0107          ! atomic weight (carbon)
        real, parameter :: massO = 15.9994          ! atomic weight (oxygen)
        real, parameter :: massN = 14.00674         ! atomic weight (nitrogen)
        real, parameter :: massS = 32.066           ! atomic weight (sulfur)
        real, parameter :: massP = 30.973761        ! atomic weight (phosphorus)
        real, parameter :: massB = 10.811           ! atomic weight (boron)
        real, parameter :: massK = 39.0983          ! atomic weight (potassium)
        real, parameter :: massF = 18.9984032       ! atomic weight (fluorine)
        real, parameter :: massI = 126.90447        ! atomic weight (iodine)

        call inputmass( natom, stmass )
        call inputfile( bffilename, frame, natom, refcoord )
        call calcmean( natom, frame, refcoord, stmass, center, cencoord )
        call outputfile( bffilename, natom, cencoord )
        allocate( ev(natom, 3) )
        call varcor( natom, cencoord, varm, frame, stmass )
        call eigvec( natom, varm )
!print *, ev(1:natom, 1)
!        call calcentropy(natom, ev, entropy, stmass)
!print *, entropy


contains
!Input the number of solute and masses of all atoms
subroutine inputmass(natom, stmass)
        implicit none
        integer :: i, j, natom
        integer, dimension(:), allocatable :: num
        character , dimension(:), allocatable :: eltp
        double precision, allocatable, dimension(:) :: stmass, charge, ljep, ljsi
        character SltInfo*20

        open(19, file = 'SltInfo', status = 'old')
          j = 1
          do
            read (19, *, end = 999)
            j = j + 1
          enddo
        999 natom = j - 1
        close(19)

        open(21, file = 'SltInfo', status = 'old')
          allocate( stmass(natom), num(natom), eltp(natom),&
                    charge(natom), ljep(natom), ljsi(natom) )
          do i = 1, natom
             read(21, *), num(i), eltp(i), charge(i), ljep(i), ljsi(i)
             if(eltp(i) == 'H') stmass(i) = massH
             if(eltp(i) == 'C') stmass(i) = massC
             if(eltp(i) == 'O') stmass(i) = massO
             if(eltp(i) == 'N') stmass(i) = massN
             if(eltp(i) == 'S') stmass(i) = massS
             if(eltp(i) == 'P') stmass(i) = massP
             if(eltp(i) == 'B') stmass(i) = massB
             if(eltp(i) == 'K') stmass(i) = massK
             if(eltp(i) == 'F') stmass(i) = massF
             if(eltp(i) == 'I') stmass(i) = massI
          enddo
        close(21)
end subroutine inputmass

!Read inputfile data
subroutine inputfile(bffilename, frame, natom, refcoord)
        implicit none
        integer :: i, frame, slvnum, natom, bfslt_natom
        double precision, allocatable, dimension(:) :: stmass
        double precision, allocatable, dimension(:,:,:) :: refcoord
        character bffilename*50

!        print *, "The number of frame is"
         read *, frame
!        print *, "The number of water molecule is"
         read *, slvnum
!        print *, "Input reffile name is"
         read *, bffilename
        call xtcbf % init( bffilename )
        call xtcbf % read
             bfslt_natom = xtcbf % NATOMS - slvnum * 3
          if(natom /= bfslt_natom) then
             print *, "Error:Inputfile is not correct"
             stop
          endif
             allocate( refcoord(3, natom, frame) )
          do i = 1, frame
             refcoord(:,:,i) = xtcbf % pos(:, 1:natom)
             call xtcbf % read
          enddo
        call xtcbf % close
end subroutine inputfile

!Calcurate mean structure
subroutine calcmean(natom, frame, refcoord, stmass, center, cencoord)
        implicit none
        integer :: i, j, k, frame, natom
        double precision :: center(3)
        double precision, allocatable, dimension(:) :: stmass
        double precision, allocatable, dimension(:,:) :: sumcoord, meancoord, outcoord, bfcoord
        double precision, allocatable, dimension(:,:,:) :: refcoord, cencoord

        call center_mass(natom, frame, refcoord, stmass, center)
        allocate( meancoord(3, natom), cencoord(3, natom, frame), &
                  sumcoord(3, natom), outcoord(3, natom), bfcoord(3, natom) )
        do j = 1, frame
           do i = 1, natom
                cencoord(1,i,j) = refcoord(1,i,j) - center(1)
                cencoord(2,i,j) = refcoord(2,i,j) - center(2)
                cencoord(3,i,j) = refcoord(3,i,j) - center(3)
           enddo
        enddo

!収束loop*5
do k = 1, 5
        do i = 1, natom
                sumcoord(1,i) = sum(cencoord(1,i,:))
                sumcoord(2,i) = sum(cencoord(2,i,:))
                sumcoord(3,i) = sum(cencoord(3,i,:))
        enddo
                meancoord(:,:) = sumcoord(1:3,1:natom) / frame
        do j = 1, frame
                bfcoord(:,:) = cencoord(1:3,1:natom,j)
                call fit(natom, meancoord, bfcoord, stmass, outcoord)
                cencoord(:,:,j) = outcoord(1:3,1:natom)
        enddo
enddo
end subroutine calcmean

!Variance-covariance matrix
subroutine varcor( natom, cencoord, varm, frame, stmass )
        integer :: i, j, k, l, natom, frame
        double precision :: cencoord(3, natom, frame), xm(natom*3), &
                            stmass(natom), nstmass(natom*3)
        double precision, allocatable :: x(:,:), varm(:,:)  !variance-covariance matrix
        allocate ( x(natom*3, frame), varm(natom*3, natom*3) )

        !Convert 3*N*frame matrix to 3N*frame matrix
        do j = 1, frame
           x(1:natom,j) = cencoord(1,1:natom,j)
           x(natom+1:natom*2,j) = cencoord(2,1:natom,j)
           x(natom*2+1:natom*3,j) = cencoord(3,1:natom,j)
        enddo
           nstmass(1:natom)=stmass(1:natom)
           nstmass(natom+1:natom*2)=stmass(1:natom)
           nstmass(natom*2+1:natom*3)=stmass(1:natom)

        !Average
           do i = 1, natom*3
              xm(i) = 0.0D0
              do j = 1, frame
                 xm(i) = xm(i) + x(i,j)*sqrt(nstmass(i))
              enddo
              xm(i) = xm(i) / dble(frame)
!print *, xm(i)
           enddo
        !Standard deviation
           do j = 1, natom*3
              do k = 1, natom*3
                 varm(j,k) = 0.0D0
                 do l = 1, frame
                 varm(j,k) = varm(j,k)+(x(j,l)-xm(j))*(x(k,l)-xm(k))
                 enddo
                 varm(j,k) = varm(j,k) / dble(frame)
!print *, varm(j,k)
              enddo
           enddo
end subroutine varcor

!eigen value analysis
subroutine eigvec( natom, varm )
        implicit none
        integer :: i, k, info, j, lda, ldvr, lwork, lwkopt, n, natom, lwmax
        complex*16 eig
        double Precision a(natom*3, natom*3), dummy(natom*3, natom*3), &
                         work(300*natom), vr(natom*3, natom*3), wi(natom*3), wr(natom*3)
        double precision :: varm(natom*3, natom*3)
        external dgeev
      
        n = natom*3
        lda=natom*3
        ldvr=natom*3
        lwork=300*natom
        a(:,:) = varm(1:natom*3,1:natom*3)
           call dgeev('N', 'N', n, a, lda, wr, wi, dummy, natom, vr, ldvr, work, lwork, info)
           lwkopt = work(1)
           if (info==0) Then
              do j = 1, n
                 write (6, *)
                 if (wi(j)==0.0D0) then
                    write (6, 100) 'Eigenvalue(', j, ') = ', wr(j)
                 else
                    eig = dcmplx(wr(j), wi(j))
                    write (6, 110) 'Eigenvalue(', j, ') = ', eig
                 endif
              enddo
           else
             write (6, *)
             write (6, 120) 'Failure in DGEEV.  INFO = ', info
           endif
           if (lwork<lwkopt) Then
             write (6, *)
             write (6, 130) 'Optimum workspace required = ', lwkopt, &
                           'Workspace provided         = ', lwork
           endif
        100 Format (1X, A, I4, A, 1P, E11.4)
        110 Format (1X, A, I4, A, '(', 1P, E11.4, ',', 1P, E11.4, ')')
        120 Format (1X, A, I4)
        130 Format (1X, A, I5, /1X, A, I5)
end subroutine eigvec

!Calculate entropy
subroutine calcentropy( natom, ev, entropy, stmass )
        implicit none
        double precision, parameter :: k = 1.38065e-23, e = 2.718282,&
                                       h = 6.62607e-34/(2*3.14159265)
        integer :: T, natom
        double precision :: s, entropy, ev(natom, 3), stmass(natom)

!        print *, "Temperature is"
        read *, T 
           s = 1
        do i = 1, 3
           do j = 1, natom
              s = s*(1+k*T*e**2*stmass(j)*ev(j,i)/h**2)
           enddo
        enddo
print *, s
           entropy = 0.5*k*log(s)
end subroutine calcentropy

!Output to xtcfile
subroutine outputfile( bffilename, natom, cencoord )
        implicit none
        integer :: i, natom
        double precision, allocatable, dimension(:,:,:) :: cencoord
        real, allocatable, dimension(:,:) :: pos
        character bffilename*50, outfilename*50

!        print *, "Outputfile name is"
        read *, outfilename
        call xtcbf % init( bffilename )
        call xtc_out % init(outfilename, 'w')
        call xtcbf % read
        do i = 1, frame
                allocate( pos(3, natom) )
                pos(:,:) = real( cencoord(1:3, 1:natom, i) )
                call xtc_out % write(natom, xtcbf % step, xtcbf % time,&
                                     xtcbf % box, pos, xtcbf % prec)
                deallocate( pos )
        enddo
        call xtc_out % close
        call xtcbf % close
end subroutine outputfile

!Calculate center of mass
subroutine center_mass(natom, frame, refcoord, stmass, center)
        implicit none
        integer, intent(in) :: natom
        double precision, intent(in) :: refcoord(3, natom, frame), stmass(natom)
        double precision, intent(out) :: center(3)
        double precision :: summass
        integer :: i, j, frame

        summass = sum(stmass)
        if(summass == 0) write(*,*), "error: sum of mass = 0"
        do j = 1, frame
           do i = 1, 3
              center(i) = dot_product(refcoord(i,:,j), stmass(:)) / summass
           end do
        enddo
end subroutine center_mass

end program bestfit_ty
