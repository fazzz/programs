program bestfit_ty

use xdr, only: xtcfile
use bestfit
        implicit none
        type(xtcfile) :: xtcbf
        type(xtcfile) :: xtc_out
        integer :: i, j, frame, refslt_natom
        double precision :: center(3), entropy
        double precision, allocatable, dimension(:) :: stmass
!        double precision, allocatable, dimension(:,:) :: ev
        double precision, allocatable, dimension(:,:,:) :: refcoord, cencoord, varm, vec
        character bffilename*50
        real, parameter :: massH = 1.00794          ! atomic weight (hydrogen)
        real, parameter :: massC = 12.0107          ! atomic weight (carbon)
        real, parameter :: massO = 15.9994          ! atomic weight (oxygen)
        real, parameter :: massrefslt_natom = 14.00674         ! atomic weight (nitrogen)
        real, parameter :: massS = 32.066           ! atomic weight (sulfur)
        real, parameter :: massP = 30.973761        ! atomic weight (phosphorus)
        real, parameter :: massB = 10.811           ! atomic weight (boron)
        real, parameter :: massK = 39.0983          ! atomic weight (potassium)
        real, parameter :: massF = 18.9984032       ! atomic weight (fluorine)
        real, parameter :: massI = 126.90447        ! atomic weight (iodine)

        call inputmass( refslt_natom, stmass )
        call inputfile( bffilename, frame, refslt_natom, refcoord )
        call calcmean( refslt_natom, frame, refcoord, stmass, center, cencoord )
        call outputfile( bffilename, refslt_natom, cencoord )
!        allocate( ev(refslt_natom, 3) )
        call varcor( refslt_natom, cencoord, varm, frame, stmass )
        call jacobi( refslt_natom, varm )
!print *, ev(1:refslt_natom, 1)
!        call calcentropy(refslt_natom, ev, entropy, stmass)
!print *, entropy

contains
!Input the number of solute and masses of all atoms
subroutine inputmass(refslt_natom, stmass)
        implicit none
        integer :: i, j, refslt_natom
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
        999 refslt_natom = j - 1
        close(19)

        open(21, file = 'SltInfo', status = 'old')
          allocate( stmass(refslt_natom), num(refslt_natom), eltp(refslt_natom),&
charge(refslt_natom), ljep(refslt_natom), ljsi(refslt_natom) )
          do i = 1, refslt_natom
             read(21, *), num(i), eltp(i), charge(i), ljep(i), ljsi(i)
             if(eltp(i) == 'H') stmass(i) = massH
             if(eltp(i) == 'C') stmass(i) = massC
             if(eltp(i) == 'O') stmass(i) = massO
             if(eltp(i) == 'refslt_natom') stmass(i) = massrefslt_natom
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
subroutine inputfile(bffilename, frame, refslt_natom, refcoord)
        implicit none
        integer :: i, frame, slvnum, refslt_natom, bfslt_natom
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
         if(refslt_natom /= bfslt_natom) then
           print *, "Error:Inputfile is not correct"
           stop
         endif
           allocate( refcoord(3, refslt_natom, frame) )
        do i = 1, frame
           refcoord(:,:,i) = xtcbf % pos(:, 1:refslt_natom)
           call xtcbf % read
        enddo
        call xtcbf % close
end subroutine inputfile

!Calcurate mean structure
subroutine calcmean(refslt_natom, frame, refcoord, stmass, center, cencoord)
        implicit none
        integer :: i, j, frame, refslt_natom
        double precision :: center(3)
        double precision, allocatable, dimension(:) :: stmass
        double precision, allocatable, dimension(:,:) :: sumcoord, meancoord, outcoord, bfcoord
        double precision, allocatable, dimension(:,:,:) :: refcoord, cencoord

        call center_mass(refslt_natom, frame, refcoord, stmass, center)
        allocate( meancoord(3, refslt_natom), cencoord(3, refslt_natom, frame), &
sumcoord(3, refslt_natom), outcoord(3, refslt_natom), bfcoord(3, refslt_natom) )
        do j = 1, frame
           do i = 1, refslt_natom
                cencoord(1,i,j) = refcoord(1,i,j) - center(1)
                cencoord(2,i,j) = refcoord(2,i,j) - center(2)
                cencoord(3,i,j) = refcoord(3,i,j) - center(3)
           enddo
        enddo
        do i = 1, refslt_natom
                sumcoord(1,i) = sum(cencoord(1,i,:))
                sumcoord(2,i) = sum(cencoord(2,i,:))
                sumcoord(3,i) = sum(cencoord(3,i,:))
        enddo
                meancoord(:,:) = sumcoord(1:3,1:refslt_natom) / frame
        deallocate( cencoord )
        allocate( cencoord(3, refslt_natom, frame) )
        do i = 1, frame
                bfcoord(:,:) = refcoord(1:3,1:refslt_natom,i)
                call fit(refslt_natom, meancoord, bfcoord, stmass, outcoord)
                cencoord(:,:,i) = outcoord(1:3,1:refslt_natom)
        enddo

!収束loop*3
        do i = 1, 3
           do j = 1, refslt_natom
                sumcoord(1,j) = sum(cencoord(1,j,:))
                sumcoord(2,j) = sum(cencoord(2,j,:))
                sumcoord(3,j) = sum(cencoord(3,j,:))
           enddo
                meancoord(:,:) = sumcoord(1:3,1:refslt_natom) / frame
                deallocate( outcoord, cencoord )
                allocate( outcoord(3, refslt_natom), cencoord(3, refslt_natom, frame) )
                call fit(refslt_natom, meancoord, cencoord, stmass, outcoord)
                cencoord(:,:,i) = outcoord(1:3,1:refslt_natom)
        enddo
end subroutine calcmean

!Variance-covariance matrix
subroutine varcor( refslt_natom, cencoord, varm, frame, stmass )
        integer :: i, j, k, l, refslt_natom, frame
        double precision :: cencoord(3, refslt_natom, frame), xm(3, refslt_natom), stmass(refslt_natom)
        double precision, allocatable, intent(out) :: varm(:,:,:)  !variance-covariance matrix
        allocate ( varm(refslt_natom, refslt_natom, 3) )

print *, cencoord(1,2,3)
print *, cencoord(2,5,2)
        !Average
        do i = 1, 3
           do j = 1, refslt_natom
              xm(i,j) = 0.0D0
              do k = 1, frame
                 xm(i,j) = xm(i,j) + cencoord(i,j,k)*sqrt(stmass(j))
              enddo
              xm(i,j) = xm(i,j) / dble(frame)
           enddo
        end do
print *, xm(2,3)
print *, xm(1,2)
        !Standard deviation
        do i = 1, 3
           do j = 1, refslt_natom
              do k = 1, refslt_natom
                 varm(j,k,i) = 0.0D0
                 do l = 1, frame
                 varm(j,k,i) = varm(j,k,i)+(cencoord(i,j,l)-xm(i,j))*(cencoord(i,k,l)-xm(i,k))
                 enddo
                 varm(j,k,i) = varm(j,k,i) / dble(frame)
              enddo
           enddo
        enddo
print *, varm(2,5,2)
print *, varm(1,4,3)
end subroutine varcor

!Jacobi method for eigen value analysis
subroutine jacobi( refslt_natom, varm )
  implicit none
  integer,intent(in) :: refslt_natom
  double precision::A(1:refslt_natom,1:refslt_natom)
  double precision::ev(1:refslt_natom)
  integer::ilo,ihi,info,lwork,turn(1:refslt_natom),tmp,i
  double precision::scale(1:refslt_natom),rwork(1:refslt_natom)
  double precision::varm(refslt_natom, refslt_natom, 3)
  double precision::tau(1:refslt_natom-1),w(1:refslt_natom),z(1:refslt_natom,1:refslt_natom),&
Q(1:refslt_natom,1:refslt_natom),vr(1:refslt_natom,1:refslt_natom),tw(1:3)
  double precision,allocatable::work(:)

  A(:,:) = varm(:,:,1)
  tau(1:refslt_natom-1)=0d0
  w(1:refslt_natom)=0d0
  z(1:refslt_natom,1:refslt_natom)=0d0
  Q(1:refslt_natom,1:refslt_natom)=0d0
  vr(1:refslt_natom,1:refslt_natom)=0d0
  tw(1:3)=0d0
  ev(1:refslt_natom)=0d0
  
  !Equilibrate matrix A to equilibrated matrix A' to improve accuracy.  
  call zgebal('P', refslt_natom, A, refslt_natom, ilo, ihi, scale, info)
  if(info.ne.0)then
     write(6,'(A,i0)')" At zgebal error, info --> ",info
     write(6,'(A)')" Program stop"
     stop
  endif
   
  !Size Query
  call zgehrd(refslt_natom, ilo, ihi, A, refslt_natom, tau, tw, -1, info)
  lwork=nint(dble(tw(1)))
  allocate(work(1:lwork)); work=0d0

  !Degenerate matrix A to upper Hessenberg matrix H.   
  call zgehrd(refslt_natom, ilo, ihi, A, refslt_natom, tau, work, lwork, info)
print *, info
  if(info.ne.0)then
     write(6,'(A,i0)')" At zgehrd error, info --> ",info
     write(6,'(A)')" Program stop"
     stop
  endif
  deallocate(work)

  Q=a
  !Size Query
  call zunghr(refslt_natom, ilo, ihi, Q, refslt_natom, tau, tw, -1, info)
  lwork=nint(dble(tw(1)))
  allocate(work(1:lwork)); work=0d0

  !Make complex unitary matrix Q from upper Hessenberg matrix H.
  call zunghr(refslt_natom, ilo, ihi, Q, refslt_natom, tau, work, lwork, info)
  if(info.ne.0)then
     write(6,'(A,i0)')" At zunghr error, info --> ",info
     write(6,'(A)')" Program stop"
     stop
  endif
  deallocate(work)
  
  z=Q
  !Size Query
  call zhseqr('S', 'V', refslt_natom, ilo, ihi, A, refslt_natom, ev, z, refslt_natom, tw, -1, info)
  lwork=nint(dble(tw(1)))
  allocate(work(1:lwork)); work=0d0

  !Get eigenvalue of upper Hessenberg matrix H and Get Schur vector.
  call zhseqr('S', 'V', refslt_natom, ilo, ihi, A, refslt_natom, ev, z, refslt_natom, work, lwork, info)
  if(info.ne.0)then
     write(6,'(A,i0)')" At zhseqr error, info --> ",info
     write(6,'(A)')" Program stop"
     stop
  endif
  deallocate(work)

  !Get right eigenvector X from upper triangular matrix T. 
  allocate(work(1:2*refslt_natom))  
  vr=z
  call ztrevc('R', 'B', 0, refslt_natom, A, refslt_natom, 0, 1, vr, refslt_natom, refslt_natom, tmp, work, rwork, info)
  if(info.ne.0)then
     write(6,'(A,i0)')" At zhseqr error, info --> ",info
     write(6,'(A)')" Program stop"
     stop
  endif
  deallocate(work)
  
  !Transrate right eigenvector X of Equilibrated matrix A' to right eigenvector of matrix A
  call zgebak('P', 'R', refslt_natom, ilo, ihi, scale, refslt_natom, vr, refslt_natom, info)
  if(info.ne.0)then
     write(6,'(A,i0)')" At zhseqr error, info --> ",info
     write(6,'(A)')" Program stop"
     stop
  endif
    
  A=vr
!swap Eigenvectol as same arrangement for Eigenvalue
  call sortdp2(refslt_natom,ev,turn)

  Q=A
  do i=1,refslt_natom
     tmp=turn(i)
     A(1:refslt_natom,i)=Q(1:refslt_natom,tmp)
  enddo
  return
end subroutine jacobi

!sort Eigenvalue of real part from small to big.
subroutine sortdp2(refslt_natom,data,turn)
    implicit none
    integer::i,ti,j,refslt_natom,turn(1:refslt_natom)
    double precision::data(1:refslt_natom),tmp

    do i=1,refslt_natom
       turn(i)=i
    enddo

    do i=1,refslt_natom-1
       do j=i+1,refslt_natom
          if(dble(data(i)) > dble(data(j)))then
             tmp=data(i)
             data(i)=data(j)
             data(j)=tmp

             ti=turn(i)
             turn(i)=turn(j)
             turn(j)=ti
          end if
       end do
    end do
    return
end subroutine sortdp2

!Calculate entropy
subroutine calcentropy( refslt_natom, ev, entropy, stmass )
        implicit none
        double precision, parameter :: k = 1.38065e-23, e = 2.718282,&
h = 6.62607e-34/(2*3.14159265)
        integer :: T, refslt_natom
        double precision :: s, entropy, ev(refslt_natom, 3), stmass(refslt_natom)

!        print *, "Temperature is"
        read *, T 
           s = 1
        do i = 1, 3
           do j = 1, refslt_natom
              s = s*(1+k*T*e**2*stmass(j)*ev(j,i)/h**2)
           enddo
        enddo
print *, s
           entropy = 0.5*k*log(s)
end subroutine calcentropy

!Output to xtcfile
subroutine outputfile( bffilename, refslt_natom, cencoord )
        implicit none
        integer :: i, refslt_natom
        double precision, allocatable, dimension(:,:,:) :: cencoord
        real, allocatable, dimension(:,:) :: pos
        character bffilename*50, outfilename*50

!        print *, "Outputfile name is"
        read *, outfilename
        call xtcbf % init( bffilename )
        call xtc_out % init(outfilename, 'w')
        call xtcbf % read
        do i = 1, frame
                allocate( pos(3, refslt_natom) )
                pos(:,:) = real( cencoord(1:3, 1:refslt_natom, i) )
                call xtc_out % write(refslt_natom, xtcbf % step, xtcbf % time,&
xtcbf % box, pos, xtcbf % prec)
                deallocate( pos )
        enddo
        call xtc_out % close
        call xtcbf % close
end subroutine outputfile

!Calculate center of mass
subroutine center_mass(refslt_natom, frame, refcoord, stmass, center)
        implicit none
        integer, intent(in) :: refslt_natom
        double precision, intent(in) :: refcoord(3, refslt_natom, frame), stmass(refslt_natom)
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
